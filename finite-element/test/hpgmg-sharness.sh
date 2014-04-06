#!/bin/bash

. ./sharness.sh

# Public: Run parallel executable and compare output
#
# When the test passed, an "ok" message is printed and the number of successful
# tests is incremented. When it failed, a "not ok" message is printed and the
# number of failed tests is incremented.
#
# With --immediate, exit test immediately upon the first failed test.
#
# Usually takes four arguments:
# $1 - Test description
# $2 - Number of processes
# $3 - Executable name (found in ${HPGMG_BINDIR}/) followed by runtime options
# $4 - Expected output
#
# With five arguments, the first will be taken to be a prerequisite:
# $1 - Comma-separated list of test prerequisites. The test will be skipped if
#      not all of the given prerequisites are set. To negate a prerequisite,
#      put a "!" in front of it.
# $2 - Test description
# $3 - Number of processes
# $4 - Executable name (found in ${HPGMG_BINDIR}/) followed by runtime options
# $5 - Expected output
#
# Returns nothing.
test_expect_stdout() {
    test "$#" = 5 && { test_prereq=$1; shift; } || test_prereq=
    test "$#" = 4 || error "bug in test script: $# not 4 or 5 parameters to test_expect_stdout"

    export test_prereq
    if ! test_skip_ "$@"; then
        say >&3 "expecting success: $2 $3"
        sed '1d;$d' <<<"$4" > reference.out
        diffoutput=
        if "${MPIEXEC}" -n $2 "${HPGMG_BINDIR}/"$3 > actual.out 2>&4 &&
            diffoutput=$(git diff --exit-code --no-index reference.out actual.out); then
            test_ok_ "$1"
        else
            test_failure_ "$1 $2 $3" "${diffoutput}"
            
        fi
    fi
    echo >&3 ""
}

# Public: Run parallel executable and check for failure with error message
#
# When the test passed, an "ok" message is printed and the number of successful
# tests is incremented. When it failed, a "not ok" message is printed and the
# number of failed tests is incremented.
#
# With --immediate, exit test immediately upon the first failed test.
#
# Usually takes four arguments:
# $1 - Test description
# $2 - Number of processes
# $3 - Executable name (found in ${HPGMG_BINDIR}/) followed by runtime options
# $4 - Expected string in error message (stderr)
#
# With five arguments, the first will be taken to be a prerequisite:
# $1 - Comma-separated list of test prerequisites. The test will be skipped if
#      not all of the given prerequisites are set. To negate a prerequisite,
#      put a "!" in front of it.
# $2 - Test description
# $3 - Number of processes
# $4 - Executable name (found in ${HPGMG_BINDIR}/) followed by runtime options
# $5 - Expected string in error message
#
# Returns nothing.
test_expect_error() {
    test "$#" = 5 && { test_prereq=$1; shift; } || test_prereq=
    test "$#" = 4 || error "bug in test script: $# not 4 or 5 parameters to test_expect_stdout"

    export test_prereq
    if ! test_skip_ "$@"; then
        say >&3 "checking known breakage: $2 $3"
        expected_stderr=$(sed '1d' <<<"$4")
        # Don't check exit code because process managers do not always propagate correctly
        "${MPIEXEC}" -n $2 "${HPGMG_BINDIR}/"$3 > /dev/null 2> actual.err
        if fgrep -q "${expected_stderr}" actual.err; then
            test_ok_ "$1"
        else
            test_failure_ "$1 $2 $3" "Expecting: ${expected_stderr}$(echo && cat actual.err)"
        fi
    fi
    echo >&3 ""
}

MPIEXEC=$(awk '/MPIEXEC/{print $3}' "${PETSC_DIR}/${PETSC_ARCH}/conf/petscvariables")
