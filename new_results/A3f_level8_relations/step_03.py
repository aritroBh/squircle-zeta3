print(open('step_01.py').read()[:3000] if __import__('os').path.exists('step_01.py') else 'no file')
