from os import environ


def setBlastDbEnv(dbDir='', origionalBLASTDBenv='', restore=False):
    '''Chenge the environment variable to where target database locates,
    call with "origionalBLASTDBenv" and "restore = True" to restore previous settings.
    This is necessary for blast program to find the database.'''
    if restore:
        print('\nRestoring BLASTDB environment variable...')
        if origionalBLASTDBenv:
            environ['BLASTDB'] = origionalBLASTDBenv
        else:
            environ.pop('BLASTDB')
        print('Restored!')

    else:
        print('\nSetting up BLASTDB environment variable...')
        if "BLASTDB" in environ:
            origionalBLASTDBenv = getenv('BLASTDB')
        else:
            origionalBLASTDBenv = False

        try:
            environ['BLASTDB'] = dbDir
            print('Success, environment variable BLASTDB will be restored in the end.\n')
        except:
            print('BLASTDB environment variable set up failed.')
            input('<Enter> to exit')
            exit()
        return origionalBLASTDBenv
# setBlastDbEnv
