#!/bin/bash

# submit ipw assays run
./submit_sl_assays_ipw.sh

# submit noise run with ipw
./submit_sl_assays_noise.sh "FALSE"

# with AIPW
./submit_sl_assays_noise.sh "TRUE"
