# get COCO dataset
mkdir dataset
mkdir dataset/COCO/
cd dataset/COCO/
git clone https://github.com/pdollar/coco.git
cd ../../

mkdir dataset/COCO/images
mkdir dataset/COCO/images/mask2014
mkdir dataset/COCO/mat
mkdir dataset/COCO/json

wget http://images.cocodataset.org/annotations/annotations_trainval2017.zip
wget http://images.cocodataset.org/zips/train2017.zip
wget http://images.cocodataset.org/zips/val2017.zip
wget http://images.cocodataset.org/zips/test2017.zip

unzip annotations_trainval2017.zip -d dataset/COCO/
unzip val2017.zip -d dataset/COCO/images
unzip test2017.zip -d dataset/COCO/images
unzip train2017.zip -d dataset/COCO/images


#rm -f annotations_trainval2014.zip
#rm -f test2015.zip
#rm -f test2014.zip
#rm -f train2014.zip
#rm -f val2014.zip
