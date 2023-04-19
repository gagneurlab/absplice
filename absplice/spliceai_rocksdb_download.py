import rocksdb
import wget
import click
import tarfile

db_url = {
    'grch37': 'https://nextcloud.in.tum.de/index.php/s/dEDnN8gBEsAbtCq/download',
    'grch38': 'https://nextcloud.in.tum.de/index.php/s/SzHjFkd9Y6jp3Dk/download',
    '_test': ''
}


@click.command()
@click.option('--version', help='SpliceAI rocksdb version (currently grch37, grch38 supported)')
@click.option('--db_path', help='Path to download database')
def spliceai_rocksdb_download(version, db_path):

    if version not in db_url:
        raise(f'Version {version} is not supported.')

    print('Downloading database...')
    download_path = db_path + '_backup.tar.gz'
    wget.download(db_url[version], out=download_path)

    print('\nUnzipping database...')
    file = tarfile.open(download_path, "r:gz")
    backup_dir = db_path + '_backup/'
    file.extractall(backup_dir)
    file.close()

    print('Storing database...')
    backup = rocksdb.BackupEngine(backup_dir + 'backup/')
    backup.restore_latest_backup(db_path, db_path)