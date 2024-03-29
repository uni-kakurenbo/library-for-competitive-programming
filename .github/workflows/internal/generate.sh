#! /bin/bash
set -eu

PROBLEM="$1"
HASH="$2"
PID="$$"

cd ./testcases/

mkdir -p "./${HASH}/"

set +e
{
    echo "::group::${PROBLEM} [HASH: ${HASH}]"

    oj download --system --silent --yukicoder-token "${YUKICODER_TOKEN}" \
        --directory "./${HASH}/" "${PROBLEM}" ||
        rm -r "./${HASH}/"

    find ~/ -type d -name "$(basename "${PROBLEM}")"
    CACHE_DIRECTORY="$(find ~/ -type d -name "$(basename "${PROBLEM}")")"

    if [ -d "${CACHE_DIRECTORY}" ]; then
        CHECKER="$(find "${CACHE_DIRECTORY}" -type f -name checker)"
        cp "${CHECKER}" "./${HASH}/checker"
    fi

    echo "::notice title=${PROBLEM}::Fetched and generated."
    echo '::endgroup::'

} &>"../log-${PID}.txt"
set -e

cat "../log-${PID}.txt"
