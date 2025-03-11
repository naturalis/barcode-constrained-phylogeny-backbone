#!/usr/bin/env python
"""
OpenToL Client Module - Interface to Open Tree of Life APIs

This module provides a client for interacting with Open Tree of Life APIs
for taxon name resolution and induced subtree retrieval.
"""

import os
import json
import time
import logging
import requests
from collections import defaultdict


class OpenToLClient:
    """Client for interacting with Open Tree of Life APIs."""

    # API base URLs
    OPENTOL_API_BASE = "https://api.opentreeoflife.org/v3"
    TNRS_MATCH_NAMES_URL = f"{OPENTOL_API_BASE}/tnrs/match_names"
    INDUCED_SUBTREE_URL = f"{OPENTOL_API_BASE}/tree_of_life/induced_subtree"
    SUBTREE_URL = f"{OPENTOL_API_BASE}/tree_of_life/subtree"
    MRCA_URL = f"{OPENTOL_API_BASE}/tree_of_life/mrca"

    # Maximum number of taxon names per batch for TNRS
    MAX_NAMES_PER_BATCH = 1000

    # Maximum number of OTT IDs per induced subtree request
    MAX_IDS_PER_SUBTREE = 500

    def __init__(self, config=None):
        """
        Initialize with optional configuration for API settings.

        Args:
            config (dict, optional): Configuration for API settings.
                                     May include 'cache_dir', 'batch_size',
                                     'rate_limit', etc.
        """
        self.config = config or {}
        self.logger = logging.getLogger(__name__)

        # Setup cache
        self.cache_dir = self.config.get('cache_dir', '.opentol_cache')
        self._setup_cache()

        # Request settings
        self.batch_size = min(
            self.config.get('batch_size', self.MAX_NAMES_PER_BATCH),
            self.MAX_NAMES_PER_BATCH
        )
        self.rate_limit = self.config.get('rate_limit', 1.0)  # seconds between requests
        self.timeout = self.config.get('timeout', 60)  # seconds
        self.retries = self.config.get('retries', 3)

        # Flag to skip OpenToL queries (for testing or offline mode)
        self.skip_opentol = self.config.get('skip_opentol', False)

        # Initialize cache
        self.cache = {
            'tnrs': {},
            'subtree': {},
            'mrca': {}
        }

        self.logger.info(f"OpenToL client initialized with cache_dir={self.cache_dir}")
        if self.skip_opentol:
            self.logger.warning("OpenToL API queries are disabled (skip_opentol=True)")

    def resolve_names(self, taxon_names, context=None):
        """
        Resolve taxon names to OTT IDs using TNRS API.

        Args:
            taxon_names (list): List of taxon names to resolve.
            context (str, optional): Taxonomic context to limit search scope.
                                     See OpenToL docs for valid contexts.

        Returns:
            dict: Mapping of taxon names to OTT IDs and metadata.
        """
        if self.skip_opentol:
            self.logger.warning("Skipping OpenToL name resolution (API queries disabled)")
            return {name: None for name in taxon_names}

        self.logger.debug(f"Resolving {len(taxon_names)} taxon names")

        # Remove duplicates while preserving order
        unique_names = []
        seen = set()
        for name in taxon_names:
            if name not in seen:
                unique_names.append(name)
                seen.add(name)

        results = {}

        # Process names in batches to avoid API limits
        for i in range(0, len(unique_names), self.batch_size):
            batch = unique_names[i:i + self.batch_size]
            self.logger.debug(f"Processing batch of {len(batch)} names (batch {i // self.batch_size + 1})")

            batch_results = self._resolve_name_batch(batch, context)
            results.update(batch_results)

            # Rate limiting
            if i + self.batch_size < len(unique_names):
                time.sleep(self.rate_limit)

        # Check for unresolved names
        unresolved = [name for name in taxon_names if name in results and results[name] is None]
        if unresolved:
            percent_unresolved = len(unresolved) / len(taxon_names) * 100
            self.logger.warning(f"Could not resolve {len(unresolved)} names ({percent_unresolved:.1f}%)")

        return results

    def get_subtree(self, node_id):
        """
        Get induced subtree from OpenToL using OTT IDs.

        Args:
            node_id (str): ID of the subtree , e.g. mrcaott170987ott201497

        Returns:
            dict: Induced subtree response containing topology information.
        """
        if self.skip_opentol:
            self.logger.warning("Skipping OpenToL subtree request (API queries disabled)")
            return None

        if not node_id:
            self.logger.warning("Node node ID provided to get_subtree")
            return None

        self.logger.info(f"Fetching subtree subtended by {node_id}")

        # Check cache first
        cache_key = self._make_cache_key('subtree', node_id)
        cached_result = self._get_from_cache('subtree', cache_key)
        if cached_result:
            self.logger.debug(f"Using cached subtree result for {node_id} root")
            return cached_result

        # Prepare request
        payload = {
            'node_id': node_id
        }

        # Make request
        return self._do_request(cache_key, payload, self.SUBTREE_URL)

    def get_induced_subtree(self, ott_ids):
        """
        Get induced subtree from OpenToL using OTT IDs.

        Args:
            ott_ids (list): List of OTT IDs.

        Returns:
            dict: Induced subtree response containing topology information.
        """
        if self.skip_opentol:
            self.logger.warning("Skipping OpenToL induced subtree request (API queries disabled)")
            return None

        if not ott_ids:
            self.logger.warning("Empty list of OTT IDs provided to get_induced_subtree")
            return None

        self.logger.info(f"Fetching induced subtree for {len(ott_ids)} OTT IDs")
        self.logger.info(f"OTT IDs: {ott_ids}")

        # If too many IDs, split into batches and combine results
        if len(ott_ids) > self.MAX_IDS_PER_SUBTREE:
            self.logger.warning(
                f"More than {self.MAX_IDS_PER_SUBTREE} OTT IDs provided to get_induced_subtree. "
                f"Splitting into multiple requests."
            )
            # This is complex and would require combining multiple subtrees
            # For now, we'll just use the first batch
            ott_ids = ott_ids[:self.MAX_IDS_PER_SUBTREE]
            self.logger.warning(f"Using only the first {self.MAX_IDS_PER_SUBTREE} OTT IDs")

        # Check cache first
        cache_key = self._make_cache_key('subtree', ott_ids)
        cached_result = self._get_from_cache('subtree', cache_key)
        if cached_result:
            self.logger.debug(f"Using cached induced subtree result for {len(ott_ids)} OTT IDs")
            return cached_result

        # Prepare request
        payload = {
            'ott_ids': ott_ids
        }

        # Make request
        return self._do_request(cache_key, payload, self.INDUCED_SUBTREE_URL)

    def _do_request(self, cache_key, payload, url):
        """
        Perform an API request with retries.
        """
        for attempt in range(self.retries):
            try:
                response = requests.post(
                    url,
                    json=payload,
                    headers={'Content-Type': 'application/json'},
                    timeout=self.timeout
                )

                if response.status_code == 200:
                    result = response.json()

                    # Cache result
                    self._save_to_cache('subtree', cache_key, result)

                    return result
                else:

                    # Failed first time, may be recoverable
                    self.logger.warning(
                        f"Failed to get subtree: {response.status_code} - {response.text}"
                    )
                    cache_key, payload, url = self._recoverable_payload(cache_key, payload, url, response)
                    if cache_key and payload and url:
                        return self._do_request(cache_key, payload, url)

            except requests.RequestException as e:
                self.logger.warning(f"Request error (attempt {attempt + 1}/{self.retries}): {str(e)}")

            # Wait before retry
            time.sleep(self.rate_limit * (2 ** attempt))
        self.logger.error(f"Failed to get subtree after {self.retries} attempts")

    def _recoverable_payload(self, cache_key, payload, url, response):

        # Handle pruned OTT IDs: If we remove them from the payload, we can retry
        # 400 - {
        #     "message": "[/v3/tree_of_life/induced_subtree] Error: node_id 'ott7851091' was not found!",
        #     "unknown": {
        #         "ott7851091": "pruned_ott_id"
        #     }
        # }
        if response.status_code == 400:
            try:
                error_data = response.json()
                if 'unknown' in error_data:

                    # Extract unknown OTT IDs and strip the 'ott' prefix
                    pruned_ids = list(error_data['unknown'].keys())
                    pruned_ids = [int(ott_id[3:]) for ott_id in pruned_ids]
                    self.logger.warning(f"Pruned OTT IDs detected, will remove from payload: {pruned_ids}")

                    # Remove pruned IDs and return the updated request
                    payload['ott_ids'] = [ott_id for ott_id in payload['ott_ids'] if ott_id not in pruned_ids]
                    return cache_key, payload, url
            except Exception as e:
                self.logger.warning(f"Failed to parse error response: {str(e)}")


    def get_mrca(self, ott_ids):
        """
        Get most recent common ancestor from OpenToL.

        Args:
            ott_ids (list): List of OTT IDs.

        Returns:
            dict: MRCA response containing node information.
        """
        if self.skip_opentol:
            self.logger.warning("Skipping OpenToL MRCA request (API queries disabled)")
            return None

        if not ott_ids or len(ott_ids) < 2:
            self.logger.warning("Need at least 2 OTT IDs to find MRCA")
            return None

        self.logger.info(f"Finding MRCA for {len(ott_ids)} OTT IDs")

        # Check cache first
        cache_key = self._make_cache_key('mrca', ott_ids)
        cached_result = self._get_from_cache('mrca', cache_key)
        if cached_result:
            self.logger.debug(f"Using cached MRCA result for {len(ott_ids)} OTT IDs")
            return cached_result

        # Prepare request
        payload = {
            'ott_ids': ott_ids
        }

        # Make request
        for attempt in range(self.retries):
            try:
                response = requests.post(
                    self.MRCA_URL,
                    json=payload,
                    headers={'Content-Type': 'application/json'},
                    timeout=self.timeout
                )

                if response.status_code == 200:
                    result = response.json()

                    # Cache result
                    self._save_to_cache('mrca', cache_key, result)

                    return result
                else:
                    self.logger.warning(
                        f"Failed to get MRCA: {response.status_code} - {response.text}"
                    )

            except requests.RequestException as e:
                self.logger.warning(f"Request error (attempt {attempt + 1}/{self.retries}): {str(e)}")

            # Wait before retry
            time.sleep(self.rate_limit * (2 ** attempt))

        self.logger.error(f"Failed to get MRCA after {self.retries} attempts")
        return None

    def clear_cache(self):
        """Clear the API response cache."""
        self.logger.info("Clearing OpenToL API cache")
        self.cache = {
            'tnrs': {},
            'subtree': {},
            'mrca': {}
        }

        # Clear cache files
        for cache_type in ['tnrs', 'subtree', 'mrca']:
            cache_dir = os.path.join(self.cache_dir, cache_type)
            if os.path.exists(cache_dir):
                for filename in os.listdir(cache_dir):
                    file_path = os.path.join(cache_dir, filename)
                    try:
                        if os.path.isfile(file_path):
                            os.unlink(file_path)
                    except Exception as e:
                        self.logger.warning(f"Failed to delete cache file {file_path}: {str(e)}")

    def _resolve_name_batch(self, names, context=None):
        """
        Resolve a batch of taxon names to OTT IDs.

        Args:
            names (list): Batch of taxon names to resolve.
            context (str, optional): Taxonomic context to limit search scope.

        Returns:
            dict: Mapping of taxon names to OTT IDs and metadata.
        """
        # Check cache first
        results = {}
        names_to_resolve = []

        for name in names:
            cache_key = self._make_cache_key('tnrs', [name, context or "None"])
            cached_result = self._get_from_cache('tnrs', cache_key)

            if cached_result:
                results[name] = cached_result
            else:
                names_to_resolve.append(name)

        # If all names were cached, return results
        if not names_to_resolve:
            return results

        # Prepare request
        payload = {
            'names': names_to_resolve
        }

        if context:
            payload['context_name'] = context

        # Make request
        for attempt in range(self.retries):
            try:
                response = requests.post(
                    self.TNRS_MATCH_NAMES_URL,
                    json=payload,
                    headers={'Content-Type': 'application/json'},
                    timeout=self.timeout
                )

                if response.status_code == 200:
                    response_data = response.json()

                    # Process results
                    name_to_matches = defaultdict(list)

                    for result in response_data.get('results', []):
                        name = result.get('name', '')
                        matches = result.get('matches', [])

                        for match in matches:
                            match_score = match.get('score', 0)
                            # Only include good matches
                            if match_score >= 0.9:
                                taxon_info = match.get('taxon', {})
                                ott_id = taxon_info.get('ott_id')

                                if ott_id:
                                    match_data = {
                                        'ott_id': ott_id,
                                        'matched_name': match.get('matched_name', ''),
                                        'score': match_score,
                                        'is_synonym': match.get('is_synonym', False),
                                        'rank': taxon_info.get('rank', '')
                                    }
                                    name_to_matches[name].append(match_data)

                    # Select best match for each name and cache
                    for name in names_to_resolve:
                        matches = name_to_matches.get(name, [])

                        if matches:
                            # Sort by score descending
                            sorted_matches = sorted(matches, key=lambda m: m['score'], reverse=True)
                            best_match = sorted_matches[0]

                            # Cache result
                            cache_key = self._make_cache_key('tnrs', [name, context or "None"])
                            self._save_to_cache('tnrs', cache_key, best_match)

                            results[name] = best_match
                        else:
                            # No matches found
                            self.logger.debug(f"No matches found for taxon name: {name}")

                            # Cache negative result
                            cache_key = self._make_cache_key('tnrs', [name, context or "None"])
                            self._save_to_cache('tnrs', cache_key, None)

                            results[name] = None

                    return results
                else:
                    self.logger.warning(
                        f"Failed to resolve names: {response.status_code} - {response.text}"
                    )

            except requests.RequestException as e:
                self.logger.warning(f"Request error (attempt {attempt + 1}/{self.retries}): {str(e)}")

            # Wait before retry
            time.sleep(self.rate_limit * (2 ** attempt))

        self.logger.error(f"Failed to resolve names after {self.retries} attempts")
        # Return None for all names that couldn't be resolved
        for name in names_to_resolve:
            results[name] = None

        return results

    def _setup_cache(self):
        """Set up cache directories."""
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

        for cache_type in ['tnrs', 'subtree', 'mrca']:
            cache_type_dir = os.path.join(self.cache_dir, cache_type)
            if not os.path.exists(cache_type_dir):
                os.makedirs(cache_type_dir)

    def _make_cache_key(self, cache_type, data):
        """
        Create a cache key for the data.

        Args:
            cache_type (str): Type of cache (tnrs, subtree, mrca).
            data: Data to create a key for.

        Returns:
            str: Cache key.
        """
        if isinstance(data, list):
            # Sort to ensure consistent keys regardless of order
            data = sorted(str(item) for item in data)
            key = '_'.join(data)
        else:
            key = str(data)

        # Hash the key if it's too long
        if len(key) > 100:
            import hashlib
            key = hashlib.md5(key.encode('utf-8')).hexdigest()

        return key

    def _get_from_cache(self, cache_type, key):
        """
        Get data from cache.

        Args:
            cache_type (str): Type of cache (tnrs, subtree, mrca).
            key (str): Cache key.

        Returns:
            object or None: Cached data if found, None otherwise.
        """
        # Check in-memory cache first
        if key in self.cache[cache_type]:
            return self.cache[cache_type][key]

        # Check on-disk cache
        cache_file = os.path.join(self.cache_dir, cache_type, key + '.json')
        if os.path.exists(cache_file):
            try:
                with open(cache_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)

                # Store in memory cache
                self.cache[cache_type][key] = data

                return data
            except Exception as e:
                self.logger.warning(f"Failed to read cache file {cache_file}: {str(e)}")

        return None

    def _save_to_cache(self, cache_type, key, data):
        """
        Save data to cache.

        Args:
            cache_type (str): Type of cache (tnrs, subtree, mrca).
            key (str): Cache key.
            data: Data to cache.
        """
        # Save to in-memory cache
        self.cache[cache_type][key] = data

        # Save to on-disk cache
        cache_file = os.path.join(self.cache_dir, cache_type, key + '.json')
        try:
            with open(cache_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
        except Exception as e:
            self.logger.warning(f"Failed to write cache file {cache_file}: {str(e)}")