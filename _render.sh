#!/bin/bash

export QUARTO_PROFILE=$1
quarto run _move-files.r
quarto render

# done in quarto.yml already 
#quarto run _remove-files.r