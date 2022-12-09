#!/bin/bash
echo "start openpose"
path=$(cat photo_path.txt)
cp $path E:/Graphics/bachelor_thesis/Code/openpose-master/examples/media/
cd E:/Graphics/bachelor_thesis/Code/openpose-master
./build/x64/Release/OpenPoseDemo.exe --image_dir examples/media/ --write_json examples/output/
echo "end openpose"
