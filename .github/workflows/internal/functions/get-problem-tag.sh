#! /bin/bash

function get_problem_tag() {
    read -r args

    URL="$args"
    TAG=""

    if [[ $URL == *judge.u-aizu.ac.jp* ]]; then
        TAG="aizu-online-judge/${URL//*id=/}"
    fi

    if [[ $URL == *judge.yosupo.jp* ]]; then
        TAG="yosupo-judge/${URL//*problem\//}"
    fi

    if [[ $URL == *yukicoder.me* ]]; then
        TAG="yukicoder/${URL//*problems\/no\//}"
    fi

    echo "${TAG,,}"
}

export -f get_problem_tag
