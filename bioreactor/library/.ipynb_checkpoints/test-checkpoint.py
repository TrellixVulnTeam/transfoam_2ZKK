import cobra
import pandas as pd
import cobra_fba as cfba
model = cfba.FBA({'model_path': '../../models/json/iJO1366.json'})

print('complete model: ', model.optimize())
with model:
    model.model.genes.b2296.knock_out()
    print('ackA knocked out: ', model.optimize())