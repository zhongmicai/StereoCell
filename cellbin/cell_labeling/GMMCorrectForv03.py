from warnings import filterwarnings
filterwarnings('ignore')
import tifffile as tifi
import cv2
import os
import math
import time
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from sklearn.mixture import GaussianMixture
from multiprocessing import Process
import gzip


def parse_head(gem):
    """
    %prog <stereomics-seq data>
    return number of header lines
    """
    if gem.endswith('.gz'):
        f = gzip.open(gem, 'rb')
    else:
        f = open(gem, 'rb')

    header = ''
    num_of_header_lines = 0
    eoh = 0
    for i, l in enumerate(f):
        l = l.decode("utf-8") # read in as binary, decode first
        if l.startswith('#'): # header lines always start with '#'
            header += l
            num_of_header_lines += 1
            eoh = f.tell() # get end-of-header position
        else:
            break
    # find start of expression matrix
    f.seek(eoh)

    return num_of_header_lines


def creat_cell_gxp(maskFile,geneFile,outpath='./', transposition=False, fileName='cellbin_gmm.txt'):
    """
    %prog <CellMask><Gene expression matrix> <output Path>

    return gene expression matrix under each cell
    """

    os.makedirs(outpath, exist_ok=True)

    print("Loading mask file...")
    mask = tifi.imread(maskFile)

    if transposition:
        mask = mask.T

    _, maskImg = cv2.connectedComponents(mask, connectivity=4)

    print("Reading data..")
    typeColumn = {"geneID": 'str',
                  "x": np.uint32,
                  "y": np.uint32,
                  "values": np.uint32,
                  "UMICount": np.uint32,
                  "MIDCount": np.uint32,
                  "MIDCounts": np.uint32}

    header = parse_head(geneFile)
    genedf = pd.read_csv(geneFile, header=header, sep='\t', dtype=typeColumn)
    if "UMICount" in genedf.columns:
        genedf = genedf.rename(columns={'UMICount':'MIDCount'})
    if "MIDCounts" in genedf.columns:
        genedf = genedf.rename(columns={'MIDCounts':'MIDCount'})

    tissuedf = pd.DataFrame()
    dst = np.nonzero(maskImg)

    print("Dumping results...")
    tissuedf['x'] = dst[1] + genedf['x'].min()
    tissuedf['y'] = dst[0] + genedf['y'].min()
    tissuedf['label'] = maskImg[dst]

    res = pd.merge(genedf, tissuedf, on=['x', 'y'], how='left').fillna(0) # keep background data
    res.to_csv(os.path.join(outpath, fileName), sep='\t', index=False)
    return res


class CellCorrection(object):
    """ Cell labeling """

    def __init__(self, mask_file, gem_file, out_path, threshold, process):
        """Initialize Cell Labeling.

        Args:
          mask_file: Binarized cell segmentation image file.
          gem_file: Gene matrix file.
          out_path: Output path of gene matrix.
          threshold: local area radius.
          process: Resource scheduling parameters, the number of processes
                that need to be specified.
        """
        self.mask_file = mask_file
        self.gem_file = gem_file
        self.out_path = out_path
        self.threshold = threshold
        self.process = process
        self.radius = 50

    def __creat_gxp_data(self, ):
        """load the gene matrix into memory, returns single cell data."""
        data = creat_cell_gxp(self.mask_file, self. gem_file, outpath=self.out_path,
                              transposition=False, fileName='nuclei_mask_profile.txt')

        if 'MIDCounts' in data.columns:
            data = data.rename(columns={'MIDCounts': 'UMICount'})
        if 'MIDCount' in data.columns:
            data = data.rename(columns={'MIDCount': 'UMICount'})

        assert 'UMICount' in data.columns
        assert 'x' in data.columns
        assert 'y' in data.columns
        assert 'geneID' in data.columns

        cell_data = data[data.label != 0].copy()
        cell_coor = cell_data.groupby('label').mean()[['x', 'y']].reset_index()

        return data, cell_data, cell_coor

    def _GMM_score(self, data, cell_coor):
        radius = self.radius

        def GMM_func(x, p_num):
            """ Single cell online GMM learning """
            t0 = time.time()
            p_data = []
            if not os.path.exists(os.path.join(self.out_path, 'bg_adjust_label')):
                os.mkdir(os.path.join(self.out_path, 'bg_adjust_label'))
            for idx, i in enumerate(x):
                if (idx + 1) % 10 == 0:
                    t1 = time.time()
                    print("proc {}: {}/{} done, {:.2f}s.".format(p_num, idx, len(x), t1 - t0))
                try:
                    clf = GaussianMixture(n_components=3, covariance_type='spherical')
                    # Gaussian Mixture Model GPU version
                    cell_test = data[(data.x < cell_coor.loc[i].x + radius) & (data.x > cell_coor.loc[i].x - radius) & (
                                data.y > cell_coor.loc[i].y - radius) & (data.y < cell_coor.loc[i].y + radius)]
                    # fit GaussianMixture Model
                    clf.fit(cell_test[cell_test.label == cell_coor.loc[i].label][['x', 'y', 'UMICount']].values)
                    cell_test_bg_ori = cell_test[cell_test.label == 0]
                    bg_group = cell_test_bg_ori.groupby(['x', 'y']).agg(UMI_max=('UMICount', 'max')).reset_index()
                    cell_test_bg = pd.merge(cell_test_bg_ori, bg_group, on=['x', 'y'])		
                    # threshold 20
                    score = pd.Series(-clf.score_samples(cell_test_bg[['x', 'y', 'UMI_max']].values))
                    cell_test_bg['score'] = score
                    threshold = self.threshold
                    cell_test_bg['label'] = np.where(score < threshold, cell_coor.loc[i].label, 0)
                     # used multiprocessing have to save result to file
                    p_data.append(cell_test_bg)
                except Exception as e:
                    print(e)
                    with open(os.path.join(self.out_path, 'error_log.txt'), 'a+') as f:
                        f.write('Cell ID: {}\n'.format(cell_coor.loc[i].label))

            out = pd.concat(p_data)
            out.drop('UMI_max', axis=1, inplace=True)
            out.to_csv(os.path.join(self.out_path, 'bg_adjust_label', '{}.txt'.format(p_num)), sep='\t', index=False)

        processes = []
        qs = math.ceil(len(cell_coor.index) / int(self.process))
        for i in tqdm(range(self.process)):

            idx = np.arange(i * qs, min((i + 1) * qs, len(cell_coor.index)))
            if len(idx) == 0: continue
            p = Process(target=GMM_func, args=(idx, i))
            p.start()
            processes.append(p)

        [pi.join() for pi in processes]

        return None

    def _GMM_correction(self, cell_data):

        bg_data = []
        error = []

        file_list = os.listdir(os.path.join(self.out_path, 'bg_adjust_label'))
        for i in tqdm(file_list):
            try:
                tmp = pd.read_csv(os.path.join(self.out_path, 'bg_adjust_label', i), sep='\t')
                bg_data.append(tmp[tmp.label != 0])
            except:
                error.append(i)
        adjust_data = pd.concat(bg_data).sort_values('score')
        adjust_data = adjust_data.drop_duplicates(subset=['geneID', 'x', 'y', 'UMICount'], keep='first').rename(columns={'score':'tag'})
        adjust_data['tag'] = '1'   #adjust
        cell_data['tag'] = '0'  #raw
        correct_data = pd.concat([adjust_data, cell_data])
        correct_data.to_csv(os.path.join(self.out_path, 'cell_mask_profile.txt'), sep='\t', index=False)

    def cell_correct(self, ):
        """ Control program of Cell labeling """
        if not os.path.exists(os.path.join(self.out_path, 'bg_adjust_label')):
            os.mkdir(os.path.join(self.out_path, 'bg_adjust_label'))

        t0 = time.time()
        data, cell_data, cell_coor = self.__creat_gxp_data()
        t1 = time.time()
        print('Load data :', (t1 - t0))
        self._GMM_score(data, cell_coor)
        t2 = time.time()
        print('Calc score :', (t2 - t1))
        self._GMM_correction(cell_data)
        t3 = time.time()
        print('Correct :', (t3 - t2))
        print('Total :', (t3 - t0))


def args_parse():
    usage = """ Usage: %s Cell expression file (with background) path, multi-process """
    arg = argparse.ArgumentParser(usage=usage)
    arg.add_argument('-m', '--mask_path', help='cell mask')
    arg.add_argument('-g', '--gem_path', help='gem file')
    arg.add_argument('-o', '--out_path', help='output path', default='./')
    arg.add_argument('-p', '--process', help='n process', type=int, default=10)
    arg.add_argument('-t', '--threshold', help='threshold', type=int, default=20)

    return arg.parse_args()


def main():
    args = args_parse()
    correction = CellCorrection(args.mask_path, args.gem_path, args.out_path, args.threshold, args.process)
    correction.cell_correct()
    

if __name__ == '__main__':
    main()

