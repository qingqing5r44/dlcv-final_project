#!/bin/bash
echo "start smplifyx"
path=$(cat photo_path.txt)
cp $path E:/Graphics/bachelor_thesis/Code/smplify-x-master/data/images
array=(${path//// })
name_part=${array[-1]%.jpg*}
name="E:/Graphics/bachelor_thesis/Code/openpose-master/examples/output/${name_part}_keypoints.json"
cp $name E:/Graphics/bachelor_thesis/Code/smplify-x-master/data/keypoints
cd E:/Graphics/bachelor_thesis/Code/smplify-x-master
source activate python37
python smplifyx/main.py --config cfg_files/fit_smpl.yaml  --data_folder data --output_folder output_smpl --visualize="False" --model_folder models --vposer_ckpt vposer --part_segm_fn smplx_parts_segm.pkl
deactivate python37
echo "end smplifyx"