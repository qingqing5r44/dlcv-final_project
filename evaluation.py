import unittest
import torch
from collections import OrderedDict
import sys  
sys.path.append(r'/home/jupyter/pytorch_Realtime_Multi-Person_Pose_Estimation')  
sys.path.append(r'/home/jupyter/OpenPose-Pytorch')
sys.path.append(r'/home/jupyter/Realtime_Multi-Person_Pose_Estimation/training/dataset/COCO')
sys.path.append(r'/home/jupyter/openpose-self')

from evaluate.coco_eval import run_eval
from lib.network.rtpose_vgg import get_model, use_vgg
from lib.network.openpose import OpenPose_Model, use_vgg
from models.rtpose_vgg import get_model
from torch import load
from models.cmu_model_torch import CmuModel

with torch.autograd.no_grad():
    # this path is with respect to the root of the project
    
    # eval our model
    model_path = '/home/jupyter/OpenPose-Pytorch/models/half_pose.pth'  # full_pose.pth
    model = get_model(trunk='vgg19')
    model = torch.nn.DataParallel(model)
    if torch.cuda.is_available():
        model = model.cuda()
    model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
    
    ## eval original model
    #model_path = '/home/jupyter/OpenPose-Pytorch/models/body_full.pth'
    #model = CmuModel()
    #if torch.cuda.is_available():
    #    model = model.cuda()
    #model.load_state_dict(torch.load(model_path))

    model.eval()
    model.float()
    model = model.cuda()
    run_eval(image_dir= '/home/jupyter/Realtime_Multi-Person_Pose_Estimation/training/dataset/COCO/images/val2017', anno_file = '/home/jupyter/Realtime_Multi-Person_Pose_Estimation/training/dataset/COCO/annotations/person_keypoints_val2017.json', vis_dir = '/home/jupyter/Realtime_Multi-Person_Pose_Estimation/training/dataset/COCO/images/val2017', model=model, preprocess='vgg')


