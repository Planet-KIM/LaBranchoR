import os
import json

def parser_config(rule):
    """
    Parameters
    ----------
    rule: str
        example: model.mercer.rebuild
    
    Returns
    -------
    target_path: str or dict
        datas folder...
    exception_dic: dic
        error message
    
    Examples
    --------
    model_mercer = parser_config(rule='model.mercer.rebuild')
    """
    try:
        config_file = os.path.join(os.path.dirname(os.path.realpath(__name__)), 'config/config.json')
        if not os.path.isfile(config_file):
            raise Exception('File Not Found')

        with open(config_file) as file:
            json_file = json.load(file)
            rule_array = rule.split('.')
            for i in range(len(rule_array)):
                rule_str = rule_array[i]
                if rule_str not in json_file.keys():
                    raise Exception('Key Not Found Error')
                json_file = json_file[rule_str]
            return json_file
    except Exception as e:
        return {'Error'  : e.args[0]}
