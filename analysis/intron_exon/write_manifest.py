"""Record SHA256 checksums for an intentionally finalized release directory."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def release_files():
    for path in sorted(ROOT.rglob('*')):
        relative = path.relative_to(ROOT)
        if not path.is_file() or path.name == 'manifest.json':
            continue
        if any(part in ('__pycache__', '.git', '.validation', '.venv',
                        '.ipynb_checkpoints') for part in relative.parts):
            continue
        if 'empirical/derived/' in relative.as_posix() or path.suffix in ('.pyc', '.zip'):
            continue
        yield path


def main():
    files = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in release_files()}
    manifest = {'format': 'SHA256',
        'purpose': 'Checksums of finalized S1 Code files, excluding this manifest and generated caches.',
        'sha256': files}
    (ROOT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Recorded {len(files)} release-file hashes.')


if __name__ == '__main__':
    main()
