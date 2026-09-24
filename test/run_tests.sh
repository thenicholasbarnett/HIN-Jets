#!/usr/bin/env bash
# compile and run every test/Test*.cxx against include/, from the repo root
# usage: test/run_tests.sh [TestName ...]

set -u

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root" || exit 1

if ! command -v root-config >/dev/null 2>&1; then
	echo "root-config not found, set up ROOT first (e.g. cmsenv)"
	exit 1
fi

cxx="$(root-config --cxx)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
export TEST_TMPDIR="$tmp"

if [[ $# -gt 0 ]]; then
	tests=()
	for name in "$@"; do
		tests+=("test/${name%.cxx}.cxx")
	done
else
	tests=(test/Test*.cxx)
fi

failed=()
for src in "${tests[@]}"; do
	name="$(basename "$src" .cxx)"
	echo "== $name"
	# shellcheck disable=SC2046
	if ! "$cxx" -O2 -Wall -Iinclude $(root-config --cflags) "$src" \
		$(root-config --libs) -o "$tmp/$name"; then
		echo "   compile failed"
		failed+=("$name")
		continue
	fi
	if ! "$tmp/$name"; then
		failed+=("$name")
	fi
done

echo
if [[ ${#failed[@]} -eq 0 ]]; then
	echo "all ${#tests[@]} passed"
	exit 0
fi
echo "failed: ${failed[*]}"
exit 1
