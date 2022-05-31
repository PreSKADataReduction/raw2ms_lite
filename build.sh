#!/usr/bin/env bash

git rev-parse HEAD >current_version
sudo singularity build casa_`git rev-parse --short HEAD`.simg raw2ms.def
