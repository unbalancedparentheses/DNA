"""Tests for gene-gene interaction (epistasis) analysis."""

from genetic_health.epistasis import evaluate_epistasis, EPISTASIS_MODELS


class TestEpistasisModels:
    def test_all_models_have_required_fields(self):
        for model_id, model in EPISTASIS_MODELS.items():
            assert 'name' in model, f"{model_id} missing name"
            assert 'genes' in model, f"{model_id} missing genes"
            assert 'conditions' in model, f"{model_id} missing conditions"
            assert isinstance(model['genes'], set)
            assert len(model['conditions']) > 0

    def test_all_conditions_have_required_fields(self):
        for model_id, model in EPISTASIS_MODELS.items():
            for i, cond in enumerate(model['conditions']):
                assert 'required' in cond, f"{model_id}[{i}] missing required"
                assert 'effect' in cond, f"{model_id}[{i}] missing effect"
                assert 'risk_level' in cond, f"{model_id}[{i}] missing risk_level"
                assert 'mechanism' in cond, f"{model_id}[{i}] missing mechanism"
                assert 'actions' in cond, f"{model_id}[{i}] missing actions"
                assert cond['risk_level'] in ('low', 'moderate', 'high')
                assert len(cond['actions']) > 0

    def test_required_genes_subset_of_model_genes(self):
        for model_id, model in EPISTASIS_MODELS.items():
            for cond in model['conditions']:
                for gene in cond['required']:
                    assert gene in model['genes'], (
                        f"{model_id}: required gene {gene} not in model genes"
                    )


class TestEvaluateEpistasis:
    def test_mthfr_comt_detected(self):
        findings = [
            {'gene': 'MTHFR', 'status': 'reduced', 'magnitude': 2},
            {'gene': 'COMT', 'status': 'slow', 'magnitude': 3},
        ]
        results = evaluate_epistasis(findings)
        assert len(results) >= 1
        names = [r['name'] for r in results]
        assert any('MTHFR' in n and 'COMT' in n for n in names)

    def test_mthfr_mtrr_detected(self):
        findings = [
            {'gene': 'MTHFR', 'status': 'reduced', 'magnitude': 2},
            {'gene': 'MTRR', 'status': 'reduced', 'magnitude': 2},
        ]
        results = evaluate_epistasis(findings)
        assert len(results) >= 1
        names = [r['name'] for r in results]
        assert any('MTHFR' in n and 'MTRR' in n for n in names)

    def test_no_match_returns_empty(self):
        findings = [
            {'gene': 'MTHFR', 'status': 'normal', 'magnitude': 0},
            {'gene': 'COMT', 'status': 'fast', 'magnitude': 1},
        ]
        assert evaluate_epistasis(findings) == []

    def test_single_gene_no_interaction(self):
        findings = [
            {'gene': 'MTHFR', 'status': 'reduced', 'magnitude': 2},
        ]
        mthfr_comt = [r for r in evaluate_epistasis(findings)
                       if 'COMT' in r['name']]
        assert len(mthfr_comt) == 0

    def test_empty_findings(self):
        assert evaluate_epistasis([]) == []

    def test_caffeine_cyp1a2_comt(self):
        findings = [
            {'gene': 'CYP1A2', 'status': 'slow', 'magnitude': 2},
            {'gene': 'COMT', 'status': 'slow', 'magnitude': 3},
        ]
        results = evaluate_epistasis(findings)
        caffeine = [r for r in results if 'Caffeine' in r['name']]
        assert len(caffeine) >= 1

    def test_caffeine_cyp1a2_adora2a(self):
        findings = [
            {'gene': 'CYP1A2', 'status': 'slow', 'magnitude': 2},
            {'gene': 'ADORA2A', 'status': 'anxiety_prone', 'magnitude': 2},
        ]
        results = evaluate_epistasis(findings)
        caffeine = [r for r in results if 'Caffeine' in r['name']]
        assert len(caffeine) >= 1

    def test_blood_pressure_ace_agt(self):
        findings = [
            {'gene': 'ACE', 'status': 'elevated', 'magnitude': 2},
            {'gene': 'AGT', 'status': 'elevated', 'magnitude': 2},
        ]
        results = evaluate_epistasis(findings)
        bp = [r for r in results if 'Blood Pressure' in r['name']]
        assert len(bp) >= 1

    def test_metabolic_fto_tcf7l2(self):
        findings = [
            {'gene': 'FTO', 'status': 'risk', 'magnitude': 2},
            {'gene': 'TCF7L2', 'status': 'risk', 'magnitude': 2},
        ]
        results = evaluate_epistasis(findings)
        metab = [r for r in results if 'Metabolic' in r['name']]
        assert len(metab) == 1

    def test_bdnf_comt_stress(self):
        findings = [
            {'gene': 'BDNF', 'status': 'reduced', 'magnitude': 2},
            {'gene': 'COMT', 'status': 'slow', 'magnitude': 3},
        ]
        results = evaluate_epistasis(findings)
        stress = [r for r in results if 'Stress' in r['name']]
        assert len(stress) == 1

    def test_iron_inflammation(self):
        findings = [
            {'gene': 'HFE', 'status': 'carrier', 'magnitude': 2},
            {'gene': 'IL6', 'status': 'high', 'magnitude': 2},
        ]
        results = evaluate_epistasis(findings)
        iron = [r for r in results if 'Iron' in r['name']]
        assert len(iron) == 1

    def test_result_structure(self):
        findings = [
            {'gene': 'MTHFR', 'status': 'reduced', 'magnitude': 2},
            {'gene': 'COMT', 'status': 'slow', 'magnitude': 3},
        ]
        results = evaluate_epistasis(findings)
        assert len(results) >= 1
        r = results[0]
        assert 'id' in r
        assert 'name' in r
        assert 'genes_involved' in r
        assert 'effect' in r
        assert 'risk_level' in r
        assert 'mechanism' in r
        assert 'actions' in r
        assert isinstance(r['actions'], list)
        assert isinstance(r['genes_involved'], dict)

    def test_high_risk_sorted_first(self):
        findings = [
            {'gene': 'MTHFR', 'status': 'reduced', 'magnitude': 2},
            {'gene': 'COMT', 'status': 'slow', 'magnitude': 3},
            {'gene': 'CYP1A2', 'status': 'slow', 'magnitude': 2},
        ]
        results = evaluate_epistasis(findings)
        if len(results) >= 2:
            risk_order = {'high': 0, 'moderate': 1, 'low': 2}
            for i in range(len(results) - 1):
                assert (risk_order.get(results[i]['risk_level'], 3)
                        <= risk_order.get(results[i + 1]['risk_level'], 3))

    def test_multiple_interactions_detected(self):
        findings = [
            {'gene': 'MTHFR', 'status': 'reduced', 'magnitude': 2},
            {'gene': 'COMT', 'status': 'slow', 'magnitude': 3},
            {'gene': 'MTRR', 'status': 'reduced', 'magnitude': 2},
            {'gene': 'CYP1A2', 'status': 'slow', 'magnitude': 2},
            {'gene': 'ACE', 'status': 'elevated', 'magnitude': 2},
            {'gene': 'AGT', 'status': 'elevated', 'magnitude': 2},
        ]
        results = evaluate_epistasis(findings)
        # MTHFR+COMT, MTHFR+MTRR, caffeine CYP1A2+COMT, BP ACE+AGT
        assert len(results) >= 4

    def test_genes_involved_has_matched_status(self):
        findings = [
            {'gene': 'MTHFR', 'status': 'severely_reduced', 'magnitude': 3},
            {'gene': 'COMT', 'status': 'slow', 'magnitude': 3},
        ]
        results = evaluate_epistasis(findings)
        mthfr_comt = [r for r in results
                       if 'Methylation-Catecholamine' in r['name']]
        assert len(mthfr_comt) == 1
        assert 'MTHFR' in mthfr_comt[0]['genes_involved']
        assert 'severely_reduced' in mthfr_comt[0]['genes_involved']['MTHFR']

    def test_apoe_mthfr_cardiovascular(self):
        findings = [
            {'gene': 'APOE', 'status': 'e3_e4', 'magnitude': 3},
            {'gene': 'MTHFR', 'status': 'reduced', 'magnitude': 2},
        ]
        results = evaluate_epistasis(findings)
        cv = [r for r in results if 'Cardiovascular' in r['name']]
        assert len(cv) == 1
