def next_filenumber():
    import fnmatch
    import os

    pattern = 'BRdata*.pkl'

    files = os.listdir('.')
    highest_number = -1

    for name in files:
        if fnmatch.fnmatch(name, pattern):
            try:
                current_number = int(name[6:-4])
                if current_number > highest_number:
                    highest_number = current_number
            except ValueError:
                continue

    return highest_filenumber()+1



def highest_filenumber():
    import fnmatch
    import os

    pattern = 'BRdata*.pkl'

    files = os.listdir('.')
    highest_number = -1

    for name in files:
        if fnmatch.fnmatch(name, pattern):
            try:
                current_number = int(name[6:-4])
                if current_number > highest_number:
                    highest_number = current_number
            except ValueError:
                continue

    return highest_number