from ultralytics import  YOLO
from matplotlib import pyplot as plt
from PIL import Image

model = YOLO('yolov8n-seg.yaml')  # build a new model from YAML
model = YOLO('yolov8n-seg.pt') 

project = r"C:\Users\Aayush Joshi\Desktop\Hek_293_project\Dataset"
name = "80_epochs-"

results = model.train(data=r'C:\Users\Aayush Joshi\Desktop\Hek_293_project\Dataset\train\config.yaml',
                      project=project,
                      name=name,
                      epochs=80,
                      patience=0,
                      imgsz=256)