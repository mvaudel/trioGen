#!/usr/bin/env bash

##
# This script extracts the variants in LD with the pruned results.
##

java -Xmx4G -cp bin/triogen-0.5.0-beta/triogen-0.5.0-beta.jar no.uib.triogen.scripts_marc.ExtractLd
