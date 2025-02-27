#!/bin/bash

#  Function to summarize resolved argument states
function summarize_sig_nrm() {
    local typ_sig="${1}"  # Type of signal computation
    local scl_fct="${2}"  # Scaling factor

    #  Determine normalization method message
    case "${typ_sig}" in
        raw|unadj|unadjusted)
            mth_nrm="No normalization; returning unadjusted signal:"
            mth_nrm+=" '--typ_sig ${typ_sig}'."
            ;;
        len|len_frag)
            mth_nrm="Performing fragment-length normalization:"
            mth_nrm+=" '--typ_sig ${typ_sig}'."
            ;;
        norm|normalized)
            mth_nrm="Generating normalized coverage (Dickson et al., Sci Rep"
            mth_nrm+=" 2023): '--typ_sig ${typ_sig}'."
            ;;
        *)
            #  Should not be possible to see this
            mth_nrm="Unknown normalization method: '--typ_sig ${typ_sig}'."
            ;;
    esac

    #  Determine scaling factor message
    if [[ -n "${scl_fct}" ]]; then
        src_scl="Custom multiplicative scaling factor(s):"
        src_scl+=" '--scl_fct ${scl_fct}'."
    else
        src_scl="No multiplicative scaling factor(s)."
    fi

    #  Print resolved argument states
    echo "#################################################"
    echo "## Summary of signal normalization and scaling ##"
    echo "#################################################"
    echo ""
    echo "- Normalization method: ${mth_nrm}"
    echo "- Scaling factor source: ${src_scl}"
    echo ""
    echo ""
}
