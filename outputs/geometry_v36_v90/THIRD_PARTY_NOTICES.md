# Third-party resources

The optional text-token measurement uses tiktoken 0.11.0, distributed under the MIT license; see `licenses/tiktoken-MIT.txt`. The two public encoding vocabularies are retained unchanged in `token_cache`, with URLs and SHA-256 hashes in `token_cache/SOURCES.json`. Runtime wheels and transitive dependencies are not bundled; installing `requirements-token.txt` retrieves them under their respective terms. The mathematics and certificate replay use Python standard-library arithmetic.

External research papers are cited by link, not redistributed. The log-Sobolev theorem is an explicit analytic dependency described in PROOF_UNIFORM_REFINE.md. No external model credentials or model API usage are included.
