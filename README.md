# AlignTrim

Stand alone version of ARTIC's fieldbioinfomatics align_trim.py

Currently under development !


## Changes to fieldbioinfomatics
- Use primalbedtools for bedfile parsing
- Standalone package manager (uv)
- Improved `find_primer()` algorithm (~10x performance gain)
- Refactored existing (3) tests to pytest
- Updates output SAM header with command 


## TODO 
- Add and implement `require-full-length` flag
    - Require the read to start / stop in a matching primer site
- More tests
    - New logic 
    - Old logic 
- Package 
    - Conda
    - pypi
- Refactor into field 
- Reduce memory on non-normalised files
    - Non-normalised files can stream out rather than wait to end.


## Installation  
```bash
uv sync 
uv run aligntrim --help
```
