# Third-party notices

The project-level MIT license does not replace third-party copyright or license notices.

## tiktoken and tokenizer data

Historical token-count experiments use the `cl100k_base` and `o200k_base` vocabularies distributed for [OpenAI tiktoken](https://github.com/openai/tiktoken). Copies appear in:

- `outputs/geometry_v36_v90/token_cache/`
- `outputs/geometry_strategy_v238_v317/token_cache/`

Each directory retains `SOURCES.json` with source URLs and file hashes. The upstream copyright is **Copyright (c) 2022 OpenAI, Shantanu Jain**. The complete [upstream MIT license](outputs/geometry_v36_v90/licenses/tiktoken-MIT.txt) is included unchanged; the [official license](https://github.com/openai/tiktoken/blob/main/LICENSE) is also available online.

The downloaded optional Python packages from the original local environment are not redistributed in this repository. Install optional dependencies using the relevant historical requirements files.

## Mathematical literature

Reports link to prior papers and projects and identify methods adapted in the implementation. Citing a method is not a claim of ownership of the cited work. The project license covers the project contribution; it does not relicense external papers or software reached through links.
