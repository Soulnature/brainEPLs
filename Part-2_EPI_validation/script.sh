#!/usr/bin/env bash

bedtools sort -i ./data/promoter_region.bed |bedtools closest -a ./data/enhancer_region.bed -b stdin -io -d > closest_promoter.bed