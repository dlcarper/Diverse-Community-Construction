#!/bin/bash
for file in *.tar.gz; do tar -zxf $file; done
rm *.tar.gz
