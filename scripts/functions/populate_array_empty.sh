#!/bin/bash

#  Populate array with <fill> (default: "NA") only if it is currently empty
#+   - If the named array does not exist, it is declared as an indexed array.
#+   - Emptiness is tested nounset-safely via the `${array[0]+x}` expansion.
#+   - Existing, non-empty arrays are left unchanged.
#+ 
#+ Usage:
#+   populate_array_empty <arr_nam> <target_len> [fill]
#+ 
#+ Arguments:
#+   arr_nam     <str>  Name of the target array variable (not a reference)
#+   target_len  <int>  Desired length to populate (must be integer ≥ 0)
#+   fill        <str>  Value to write into each element (default: "NA")
#+ 
#+ Returns:
#+   0 on success, non-zero on validation or runtime error
#+
#+ Notes:
#+   - This function never expands the array variable unguarded; it uses
#+     `declare -p` to probe existence and `${array[0]+x}` to test emptiness.
#+   - Preserves any pre-existing contents if the array is already non-empty.
#+
#+ Example:
#+   ```
#+   # If arr_foo is unset or empty, make it 4×"NA"
#+   populate_array_empty "arr_foo" 4
#+
#+   # If arr_bar is unset or empty, make it 3×"0"
#+   populate_array_empty "arr_bar" 3 "0"
#+   ```
# shellcheck disable=SC2034
populate_array_empty() {
    local arr_nam="${1}"   # Name of target array variable (string)
    local target="${2}"    # Desired length if empty (int >= 0)
    local fill="${3:-NA}"  # Value to insert
    local i                # Loop index
    local probe=""         # Holds `${array[0]+x}` expansion safely

    #  Check that target length is int >=0
    if ! [[ "${target}" =~ ^[0-9]+$ ]]; then
        echo "Error: Target not integer: '${target}'" >&2
        return 1
    fi

    #  Ensure the target exists as an indexed array (nounset-safe)
    if ! eval 'declare -p '"${arr_nam}"' >/dev/null 2>&1'; then
        eval "declare -a ${arr_nam}=()"
    fi

    #  Nounset-safe emptiness probe: `${array[0]+x}` expands to "x" if element
    #+ 0 exists; if not, `${array[0]+x}` expands to empty string ("")
    eval 'probe=${'"${arr_nam}"'[0]+x}'

    #  Only populate if the array is empty
    if [[ -z "${probe}" ]]; then
        for (( i=0; i<target; i++ )); do
            eval "${arr_nam}[${i}]=\${fill}"
        done
    fi
}
