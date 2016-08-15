# Download file utility

def download_file(base_url, fh, download_dir):
    """
    Check if file is already downloaded, and download if not

    Parameters:
    base_url: the base http address of file to download
    fh: the file located at the base_url to download (will be name of download)
    download_dir: the location to save fh

    Output:
    Nothing - will save base_url/fh to download_dir/fh
    """
    import os
    from urllib.request import urlretrieve

    path = os.path.join(download_dir, fh)
    if not (os.path.exists(path)):
        urlretrieve(base_url + fh, path)
