# BUILD.md

## see https://delight.readthedocs.io/en/latest/install.html

## Via conda or pip
      pip install -r requirements.txt


## Make it
   python setup.py build_ext --inplace
  python setup.py install

## Clean
python setup.py install --record files.txt
sudo rm $(cat files.txt)
