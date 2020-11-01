# Define path
import platform
if platform.system() == 'Darwin':
	fileReader = ['open', '-a', 'TextEdit']
elif platform.system() == 'Linux':
	fileReader = ['less']
else: # 'Windows'
	fileReader = ['notepad']

def clear(file_path):
	import pathlib
	if pathlib.Path(file_path).is_file():
		with open(file_path, 'w') as clear_handle:
			write(f'Previous {file_path} cleaned\n', file_path)

	
def write(string, file_path, end = None):
	with open(file_path, 'a') as write_handle:
		print(string, end = end, file = write_handle)
		print(string, end = end)

def open_in_notepad(file_path):
	import subprocess
	try:
		run_notepad = subprocess.run([*fileReader, file_path])
	except:
		print(run_notepad)
	
def clean_str(str):
	import re
	return re.sub('[^\w_.)( -]', '.', str)
