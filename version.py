import subprocess

git_dir = "./"

try:
    pkg_version = (
        subprocess.check_output(["git", "describe"], cwd=git_dir)
        .strip()
        .decode("utf-8")
    )
except:
    pkg_version = "--"
