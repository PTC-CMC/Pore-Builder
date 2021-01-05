#!/bin/sh

. /opt/conda/etc/profile.d/conda.sh
conda activate base
conda activate pore37

if [ "$@" == "none" ]; then
	bash
else
	$@
fi
