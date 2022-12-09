#!/bin/bash
path=$(cat photo_path.txt)
array=(${path//// })
name_part=${array[-1]%.jpg*}
name="E:/Graphics/bachelor_thesis/Code/smplify-x-master/output_smpl/results/${name_part}/coef.txt"
cp $name E:/Graphics/bachelor_thesis/Code/humanbody/humanbody/generatedata/16