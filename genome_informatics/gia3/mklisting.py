import os

files = os.listdir('.')
for f in files:
	if f[-1] == 'R': lang = 'R'
	elif f[-2:] == 'py': lang = 'python'
	else: print('filetype not recognized for file', f)
	print('\large{File: ' + '\_'.join(f.split('_')) + '}')
	print('\lstinputlisting[language=python]{../code/'+f+'}\n')



#\large{File: parse\_pepvar.py}
#\lstinputlisting[language=python]{../code/parse_pepvar.py}
