from get_template import template_extractor
from rxnmapper import BatchedMapper
import json
import pandas as pd
import chardet
from func_timeout import func_timeout, FunctionTimedOut
from joblib import Parallel, delayed


def _extract(reaction,radius=1):
    try:
        return template_extractor.extract_from_reaction(reaction,radius)
    except KeyboardInterrupt:
        print('Interrupted')
        raise KeyboardInterrupt
    except Exception:
        return {'reaction_id': reaction['_id']}

def extract(reaction,radius):
    try:
        return func_timeout(20, _extract, args=(reaction,radius))
    except FunctionTimedOut:
        print('Timeout')
        return {'reaction_id': reaction['_id']}

def preprocess_data(data, radius=1):
    """
    Preprocess the data by extracting templates from the reactions and removing chirality and mapping numbers
    :param data: list of reactions
    :param radius: radius for template extraction
    :return: preprocessed data
    """
    reactant = data['reactant']
    product = data['product']
    reaction = []
    for i in range(len(reactant)):
        reaction.append(reactant[i] + '>>' + product[i])
    mapper = BatchedMapper(batch_size=20)
    reaction = list(mapper.map_reactions(reaction))
    data['reaction'] = reaction

    split_smiles = data['reaction'].str.split('>>',expand=True)
    data['reactants'] = split_smiles[0]
    data['products'] = split_smiles[1]
    
    reactions = data[['_id', 'reactants', 'products']].to_dict('records')
    print('Extracting templates...')
    templates_r1 = Parallel(n_jobs=-1, verbose=4)(delayed(extract)(reaction,1) for reaction in reactions)
    templates_r0 = Parallel(n_jobs=-1, verbose=4)(delayed(extract)(reaction,0) for reaction in reactions)
    
    templates_list1 = []
    for template in templates_r1:
        try:
            templates_list1.append(template['reaction_smarts'])
        except:
            templates_list1.append('')

    data['template_r1'] = templates_list1

    templates_list0 = []
    for template in templates_r0:
        try:
            templates_list0.append(template['reaction_smarts'])
        except:
            templates_list0.append('')

    data['template_r1'] = templates_list1
    data['template_r0'] = templates_list0

    out_data = data[['_id','reactant','product','reaction','condition','template_r0','template_r1']]

    return out_data


def main(args):
    with open(args.input_path, 'rb') as f:
        enc = chardet.detect(f.read())
    data = pd.read_csv(args.input_path, encoding=enc['encoding'])
    data = preprocess_data(data, radius=1)
    data.to_csv(args.input_path, index=False)  
    template_list = [i for i in data['template_r1'] if str(i) != 'NaN' and i != '' and str(i) != 'nan']
    template_list = [i for i in template_list if '.' not in i.split('>>')[0]]
    if args.add_radius_0:
        template_list += [i for i in data['template_r0'] if str(i) != 'NaN' and i != '' and str(i) != 'nan']
        template_list = [i for i in template_list if '.' not in i.split('>>')[0]]
    with open(args.output_path, 'w') as f:
        json.dump(template_list, f)
    template_condition = {}
    for template in template_list:
        condition = data[data['template_r1'] == template]['condition'].values
        condition = [i for i in condition if str(i) != 'NaN' and i != '' and str(i) != 'nan']
        template_condition[template] = condition
    with open(args.condition_path, 'w') as f:
        json.dump(template_condition, f)



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', type=str, default='./template/preprocessed_data.csv', help='Path to the input CSV file')
    parser.add_argument('--output_path', type=str, default='./reaction_template.json', help='Path to the output json file')
    parser.add_argument('--condition_path', type=str, default='./template_condition.json', help='Path to the condition templates json file')
    parser.add_argument('--add_radius_0', action='store_true', help='Whether to add radius 0 templates to the output')
    args = parser.parse_args()
    main(args)