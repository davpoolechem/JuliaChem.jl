#!/bin/bash

export ARPA_DIR="`pwd`/example_inputs/Water/adenine_thymine_2"

function submit_job {
  rm ${ARPA_DIR}/$1.output
  rm ${ARPA_DIR}/$1.cobaltlog
  rm ${ARPA_DIR}/$1.error
  
  qsub -O ${ARPA_DIR}/$1 -A CHM135 -q skylake_8180 -t 360 -n 1 ${ARPA_DIR}/$1.sh   
}

#submit_job adenine_thymine_2_1thread
submit_job adenine_thymine_2_7thread
#submit_job adenine_thymine_2_14thread
#submit_job adenine_thymine_2_28thread
#submit_job adenine_thymine_2_42thread
#submit_job adenine_thymine_2_56thread
#submit_job adenine_thymine_2_84thread
#submit_job adenine_thymine_2_111thread
#submit_job adenine_thymine_2_112thread

#submit_job fig1a_56thread
#submit_job fig1b_56thread
#submit_job fig1c_56thread
#submit_job fig1e_56thread
#submit_job fig1g_56thread
#submit_job fig1h_56thread
#submit_job fig1i_56thread
