import utils
from models import LogisticRegression, DistributionImputator, DecisionTree, ScoringModel, NomogramModel
import data_utils as du
import model_ensemble as me
import numpy as np
import model_evaluate as eval
import logging


def load_model(json_file):
    model_data = utils.load_json_data(json_file)
    if model_data['model_type'] == 'logistic_regression':
        model = LogisticRegression(model_data)
    elif model_data['model_type'] == 'decision_tree':
        model = DecisionTree(model_data)
    elif model_data['model_type'] == 'scoring':
        model = ScoringModel(model_data)
    elif model_data['model_type'] == 'nomogram':
        model = NomogramModel(model_data)
    else:
        raise Exception('model type [{0}] not recognised'.format(model_data['model_type']))
    logging.info('{0} loaded as a {1} model'.format(json_file, model.model_type))
    return model


def test_ensemble(model, x, outcome='death', threshold=0.5):
    """
    test ensemble method
    :param model:
    :param x:
    :param outcome:
    :param threshold:
    :return:
    """
    x = x.loc[x[outcome].notna()].copy()
    y = x[outcome].to_list()
    if model.mode in [me.VoteMode.average_score, me.VoteMode.max_score]:
        probs = model.predict_probs(x)
        result = eval.evaluate_pipeline(y, probs, model_name='ensemble model', threshold=threshold)
        return result
    else:
        return None


def test_single_model(model, x, outcome=None):
    """
    test a single model
    :param model:
    :param x:
    :param outcome:
    :return:
    """
    if outcome is None:
        outcome = model.outcome
    x = x.loc[x[outcome].notna()].copy()
    dist = model.model_data['cohort_variable_distribution']
    di = DistributionImputator(dist)
    x = di.impute(x, variables=[k for k in dist])
    predicted_probs = np.array(model.predict_prob(x))
    y = x[outcome].to_list()
    result = eval.evaluate_pipeline(y, predicted_probs)
    return model, result


def test_models_and_ensemble(model_files, x, weights=None, outcome='death', threshold=0.5, result_csv=None):
    """
    do tests on individual models and also ensemble methods
    :param model_files:
    :param x:
    :param weights:
    :param outcome:
    :param threshold:
    :param result_csv:
    :return:
    """
    results = {}
    ve = me.BasicEnsembler()
    for idx in range(len(model_files)):
        mf = model_files[idx]
        m = load_model(mf)
        m, result = test_single_model(m, x, outcome=outcome)
        ve.add_model(m, 1 if weights is None else weights[idx])
        results['{0}\n({1})'.format(m.id, m.model_type)] = result

    ve.mode = me.VoteMode.average_score
    result = test_ensemble(ve, x, threshold=threshold, outcome=outcome)
    results['ensemble model'] = result
    result_df = eval.format_result(results)
    if result_csv is not None:
        result_df.to_csv(result_csv, sep='\t', index=False)


def populate_col_by_or(x, cols, new_col_name):
    cm = []
    for idx, r in x.iterrows():
        v = 0
        for c in cols:
            if r[c] == 1:
                v = 1
                break
        cm.append(v)
    x[new_col_name] = cm
    return x


def do_test(config_file):
    """
    do the tests by using configuration file
    :param config_file:
    :return:
    """
    config = utils.load_json_data(config_file)
    partial_to_saturation_col = None if 'partial_to_saturation_col' not in config \
        else config['partial_to_saturation_col']
    x = du.read_data(config['data_file'],
                     column_mapping=config['mapping'],
                     partial_to_saturation_col=partial_to_saturation_col)
    if 'comorbidity_cols' in config:
        populate_col_by_or(x, config['comorbidity_cols'], new_col_name='comorbidity')
    model_files = config['model_files']
    for outcome in config['outcomes']:
        logging.info('testing for outcome [{0}] with #{1} models'.format(outcome, len(model_files)))
        result_file = '{0}/{1}_result.tsv'.format(config['result_tsv_folder'], outcome)
        test_models_and_ensemble(model_files,
                                 x,
                                 weights=config['weights'][outcome],
                                 outcome=outcome,
                                 threshold=0.5,
                                 result_csv=result_file
                                 )
        logging.info('result saved to {0}'.format(result_file))


if __name__ == "__main__":
    utils.setup_basic_logging(log_level='INFO', file='ensemble.log')
    do_test('./test/test_config.json')