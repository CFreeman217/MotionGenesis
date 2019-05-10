import os
[os.rename(f, f.replace('Exam','freeman_c_problem')) for f in os.listdir('.') if not f.endswith('.py')]

