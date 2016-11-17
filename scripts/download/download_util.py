# Download file utility


def download_file(base_url, filename, download_dir):
    """
    Check if file is already downloaded, and download if not

    Parameters:
    base_url: the base http address of file to download
    filename: the file located at the base_url (will be name of download)
    download_dir: the location to save fh

    Output:
    Nothing - will save base_url/fh to download_dir/fh
    """
    import os
    from urllib.request import urlretrieve

    path = os.path.join(download_dir, filename)
    if not (os.path.exists(path)):
        urlretrieve(base_url + filename, path)


def process_repeats(base_dir, filename, genome):
    """
    Process the irregularly structured repeatmasker file
    """
    from subprocess import call

    full_loc = base_dir + filename
    call(['gunzip', full_loc])

    # Replace variable length space delimiters with tabs and output tsv file
    CMD = "sed 's/ \+/\\t/g' " + full_loc[:-3] + " > " + full_loc[:-3] + ".tsv"
    call(CMD, shell=True)


def untar_mouse_domain(tar_path):
    """
    Mouse domains are downloaded as a tar archive - extracts contents
    """
    import tarfile

    tar = tarfile.open(tar_path)
    tar.extractall()
