import ROOT
import xml.etree.ElementTree as ET
import pandas as pd
import uproot
import array

def get_variable_names_from_xml(xml_file):
    """
    Extract feature variable names from a TMVA XML file.

    Parameters:
        xml_file (str): Path to the TMVA XML file

    Returns:
        list: List of feature variable names
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    variables_node = root.find('Variables')
    return [var.get('Expression') for var in variables_node.findall('Variable')], [var.get('Type') for var in variables_node.findall('Variable')]

def load_tmva_model(model_file, feature_names, feature_types):
    """
    Load a TMVA model and return the reader with variables.

    Parameters:
        model_file (str): Path to the TMVA model file
        feature_names (list): List of feature variable names
        feature_types (list): List of feature variable types

    Returns:
        tuple: (reader, variables)
    """
    reader = ROOT.TMVA.Reader()
    variables = {}

    for name, var_type in zip(feature_names, feature_types):
        # Create a mutable buffer using array
        variables[name] = array.array('f' if var_type == 'F' else 'i', [0])
        reader.AddVariable(name, variables[name])

    # Load the model
    reader.BookMVA(
        'DijetBDT',
        model_file
    )

    return reader, variables

# Example usage
if __name__ == "__main__":
    # Load input data
    input_data = uproot.open('/eos/project/h/htozg-dy-privatemc/rzou/bdt/VBF_pinnacles.root')
    input_signal = input_data['TreeS'].arrays(library='pd')
    input_background = input_data['TreeB'].arrays(library='pd')

    print(input_signal)
    print(input_background)

    output_signal_all, output_background_all = pd.DataFrame(), pd.DataFrame()

    for fold in range(0, 4):
        output_signal = input_signal[input_signal['part'] == fold]
        output_background = input_background[input_background['part'] == fold]

        xml_model_file = '/eos/project/h/htozg-dy-privatemc/rzou/bdt/dataset_TWOJET_pinnacles_fix/weights/TMVAClassification_BDT_{}.weights.xml'.format(fold)

        # Extract feature variable names
        feature_names, feature_types = get_variable_names_from_xml(xml_model_file)
        print("Feature variable names:", feature_names)
        print("Feature variable types:", feature_types)

        # Load model and variables
        try:
            reader, variables = load_tmva_model(xml_model_file, feature_names, feature_types)
        except FileNotFoundError:
            print(f"Error: File {xml_model_file} not found")
            exit(1)
        except Exception as e:
            print(f"Error: {e}")
            exit(1)

        # Make predictions
        predictions = []
        for _, row in output_signal.iterrows():
            # Update variable values
            for name in feature_names:
                variables[name][0] = row[name]
            
            # Calculate prediction
            prediction = reader.EvaluateMVA('DijetBDT')
            predictions.append(prediction)

        print(predictions)
        output_signal['bdt_score'] = predictions

        predictions = []
        for _, row in output_background.iterrows():
            # Update variable values
            for name in feature_names:
                variables[name][0] = row[name]
            
            # Calculate prediction
            prediction = reader.EvaluateMVA('DijetBDT')
            predictions.append(prediction)

        # print(predictions)
        output_background['bdt_score'] = predictions

        output_signal_all = pd.concat([output_signal_all, output_signal])
        output_background_all = pd.concat([output_background_all, output_background])

    with uproot.recreate('/eos/user/j/jiehan/root/output_rui_sample/signal.root') as f:
        f['signal'] = output_signal_all.to_dict(orient='list')

    with uproot.recreate('/eos/user/j/jiehan/root/output_rui_sample/background.root') as f:
        f['background'] = output_background_all.to_dict(orient='list')
    