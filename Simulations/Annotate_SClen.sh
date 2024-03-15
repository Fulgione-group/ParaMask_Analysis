#!/bin/bash
awk '{printf"%s\t%s\n", $1,"SC_b"NR}' $1
