from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

download_base = Path('/app/data')
abc_cache = AbcProjectCache.from_cache_dir(download_base)

print(abc_cache.current_manifest)