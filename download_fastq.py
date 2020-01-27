import argparse
import glob
import hashlib
import os
import random
from urllib import request
import sys
import pandas as pd
import requests
from tqdm import tqdm


# TODO: log skipped studies in file
# TODO: add test run capability
# TODO: load .csv instead of .xlsx

# straight from tqdm example:
def progress_hook(t):
    """Wraps tqdm instance.
    Don't forget to close() or __exit__()
    the tqdm instance once you're done with it (easiest using `with` syntax).
    Example
    -------
    >>> with tqdm(...) as t:
    ...     reporthook = my_hook(t)
    ...     urllib.urlretrieve(..., reporthook=reporthook)
    """
    last_b = [0]

    def update_to(b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            t.total = tsize
        t.update((b - last_b[0]) * bsize)
        last_b[0] = b

    return update_to


def download_file(url, path):

    file_dir = os.path.dirname(path)
    if not os.path.isdir(file_dir):
        os.makedirs(file_dir)
    if os.path.isfile(path):
        print("Skipping {0} (already exists)".format(path))
    else:
        with tqdm(unit="B", unit_scale=True, miniters=1, desc=os.path.basename(path)) as t:
            reporthook = progress_hook(t)
            # print('downloading', url, path)
            # shutil.copy2(url, path)
            request.urlretrieve(url, path, reporthook=reporthook)

def download_file_metadata(study_ids, output_dir, refresh=False):
    for study_id in study_ids:
        study_dir = os.path.join(output_dir, study_id)
        study_file = os.path.join(study_dir, '{0}.tsv'.format(study_id))
        # print(study_file, os.path.isfile(study_file))
        print('{0}: checking study meta data...'.format(study_id))
        if not os.path.isfile(study_file) or refresh:
            print('{0}: downloading meta data...'.format(study_id))
            # request.urlretrieve('https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={0}&result=read_run&fields='
            metadata = requests.get('https://www.ebi.ac.uk/ena/data/warehouse/filereport?'
                                    'accession={0}&'
                                    'result=read_run&'
                                    'fields='
                                    'run_accession,'
                                    'fastq_ftp,'
                                    'fastq_md5,'
                                    'fastq_bytes,'
                                    'read_count,'
                                    'base_count,'
                                    'library_layout'.format(study_id)).text.strip()
            if not os.path.isdir(study_dir):
                os.makedirs(study_dir)

            with open(study_file, 'w') as file:
                file.write(metadata)


def get_run_accessions(study_id, output_dir):
    return pd.read_csv(os.path.join(output_dir, study_id, study_id+'.tsv'), delimiter='\t')


def check_layout(accession):
    ftp_files = accession['fastq_ftp'].split(';') if not pd.isna(accession['fastq_ftp']) else []
    layout = accession['library_layout']
    name = accession['run_accession']

    if layout == 'PAIRED' and len(ftp_files) == 2:
        filenames = [os.path.splitext(os.path.splitext(os.path.basename(filename))[0])[0] for filename in ftp_files]

        if name+'_1' in filenames and name+'_2' in filenames:
            return True
        else:
            return False
    elif layout == 'SINGLE' and len(ftp_files) == 1:
        return True
    else:
        return False


def check_fastq_md5(fastq_file, fastq_md5):
    hash_md5 = hashlib.md5()
    with open(fastq_file, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest() == fastq_md5


def delete_file(file):
    os.remove(file)
    print('Deleting {0}...'.format(file))


def download_accession(acc_id, study_files, existing_fastq, output_dir):
    print('Downloading files for run {0} to {1}'.format(acc_id, output_dir))
    files = study_files[(study_files.accession == acc_id) & (study_files.file_ok == False)]
    for i, file in files.iterrows():

        download_file('ftp://'+file['ftp'], os.path.join(output_dir, i))

        study_files.loc[i, 'file_ok'] = True
        existing_fastq.append(i)


def select_runs(study_id, accessions, output_dir, max_files):
    study_dir = os.path.join(output_dir, study_id)

    study_files = {}
    for i, accession in accessions.iterrows():
        fastq_ftp = accession['fastq_ftp'].split(';')
        fastq_files = [os.path.basename(fastq) for fastq in fastq_ftp]
        fastq_md5 = accession['fastq_md5'].split(';')

        # print(fastq_files, fastq_ftp, fastq_md5)
        for i, file in enumerate(fastq_files):
            # print(i, file)
            study_files[file] = [accession['run_accession'], fastq_ftp[i], fastq_md5[i]]
    study_files = pd.DataFrame.from_dict(study_files, orient='index', columns=['accession', 'ftp', 'md5'])
    study_files['file_ok'] = False

    existing_fastq = [fastq for fastq in glob.glob(os.path.join(study_dir, '*.fastq.gz'))]
    existing_fastq = [os.path.basename(fastq) for fastq in existing_fastq]

    print('{0}: checking local fastq files...'.format(study_id))
    for file in list(existing_fastq):

        file_path = os.path.join(output_dir, study_id, file)
        if file in study_files.index:
            if not check_fastq_md5(file_path, study_files.loc[file, 'md5']):
                existing_fastq.pop(existing_fastq.index(file))
                print('{0}: local fastq file {1} did not pass MD5 check.'
                      'The file is deleted and will be redownloaded'.format(study_id, file))
                delete_file(file_path)

            else:
                study_files.loc[file, 'file_ok'] = True
                print('{0}: {1} ok'.format(study_id, file))
        else:
            existing_fastq.pop(existing_fastq.index(file))
            print('{0}: local fastq file {1} is not in the found in run accessions metadata.'
                  'The file will be ignored, but might be used by downstream programs'.format(study_id, file))


    # check for incomplete accessions
    incomplete_accessions = []
    for i, accession in accessions.iterrows():
        acc_id = accession['run_accession']
        file_ok = study_files.loc[study_files.accession == acc_id, 'file_ok']

        # if not all files are completed but at least one is
        if not file_ok.all() and file_ok.any():
            incomplete_accessions.append(acc_id)
    if len(incomplete_accessions) > 0:
        print('{0}: incomplete accessions found locally, downloading missing files...'.format(study_id))

    for acc_id in incomplete_accessions:
        download_accession(acc_id, study_files, existing_fastq, os.path.join(output_dir, study_id))

    # print(len(existing_fastq), max_files, study_files['file_ok'].all() == False)
    if len(existing_fastq) < max_files and study_files['file_ok'].all() == False:
        print('{0}: downloading remaining files...')
        while len(existing_fastq) < max_files and study_files['file_ok'].all() == False:
            remaining_accessions = study_files.loc[study_files['file_ok'] == False, 'accession'].unique()
            download_accession(random.choice(remaining_accessions), study_files, existing_fastq, os.path.join(output_dir, study_id))
        print('{0}: finished downloading.'.format(study_id))
    print('{0}: downloaded {1} of {2} available files.'.format(study_id, sum(study_files['file_ok']), len(study_files)))
    if len(existing_fastq) >= max_files:
        print('{0}: maximum number of files reached.'.format(study_id))

def download_fastq(args_in):

    parser = argparse.ArgumentParser(description="Download fastq files from metagenomic studies in the ENA database. "
                                                 "'selected_studies.txt' in the output folder "
                                                 " is the file with study accessions. "
                                                 "Meta data on the run files (fastq) are downloaded for each study and "
                                                 "placed in the study subdirectory. The meta data is checked to "
                                                 "determine if studies are pair-ended or single-ended and studies with "
                                                 "inconsistent or mixed files are ignored. Downloading can be forced by "
                                                 " placing the study IDs in 'force_download.txt', separated by a new "
                                                 "line in the output folder."
                                                 "The files are downloaded per accession number. A maximum number of "
                                                 "files per study can be specified, and new accessions wil be "
                                                 "downloaded randomly until this maximum number is passed."
                                                 "If fastq files are already present in the output folder, they will be"
                                                 "checked against the checksum in the meta data. If the checksums do "
                                                 "not match, the local file is deleted and downloaded, else it will be "
                                                 "skipped. Accession for which some files are not present will be "
                                                 "completed regardless of max file number. As long as the "
                                                 "max file number is not exceeded, additional runs will be downloaded.")

    parser.add_argument("-o", "--output_dir", dest="output_dir", required=True,
                        metavar="[output_dir]",
                        help="destination of files to be downloaded. "
                             "This directory should also contain the file selected_studies.txt as input")
    parser.add_argument("-n", '--max_files', dest="max_files", default=10, type=int,
                        metavar="[max_files]",help="specify a maximum number of files per study")
    parser.add_argument('--refresh_metadata', action='store_true', help="force downloading of run meta data")

    args = parser.parse_args(args_in)

    output_dir = args.output_dir
    max_files = args.max_files
    refresh_metadata = args.refresh_metadata

    # Check if the output directory exists
    if not os.path.isdir(output_dir):
        raise FileNotFoundError('Output directory "{0}"  could not be found'.format(output_dir))

    # Create a list of study ids from the study catalog
    studies_file = os.path.join(output_dir, 'selected_studies.txt')

    with open(studies_file, 'r') as file:
        study_list = file.read().strip().split('\n')
    print("Selected studies:\n{0}".format('\n'.join(study_list)))

    ignore_file = os.path.join(output_dir, 'force_download.txt')
    if os.path.isfile(ignore_file):
        with open(ignore_file, 'r') as file:
            ignore_layout = file.read().split('\n')
    else:
        ignore_layout = []

    download_file_metadata(study_list, output_dir, refresh=refresh_metadata)

    for study in study_list:
        accessions = get_run_accessions(study, output_dir)
        pairs = accessions['library_layout'].unique()
        if len(pairs) == 1 or study in ignore_layout:
            accessions['layout_check'] = accessions.apply(check_layout, axis=1)
            if accessions['layout_check'].all():
                select_runs(study, accessions, output_dir, max_files)
            else:
                print('{0}: mismatch between the layout type and the expected number of files in one or more accessions. '
                      'Add study id to "force_download.txt" to include.'.format(study))
        else:
            print('{0}: both single and pair-ended library layouts or unrecognised layout. '
                  'Add study id to "force_download.txt" to include.'.format(study))


if __name__ == '__main__':
    download_fastq(sys.argv[1:])

# output_dir = r"D:\17 Dutihl Lab\_tools\get_meta"
# download_fastq(["-o", output_dir, "-n", "2"])