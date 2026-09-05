"""
Security and input validation tests for vcf-consensus-heterozygosity-rate.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
import os
import warnings
from vcf_titv import _validate_path, calculate_metrics, process_batch
from agents.base import AuditTrail, PHIGuard, SecurityException


class TestPathValidation:
    """Tests for path traversal prevention."""

    def test_normal_path_passes(self):
        assert _validate_path("data/input.csv") == "data\\input.csv" or _validate_path("data/input.csv") == "data/input.csv"

    def test_simple_filename_passes(self):
        assert _validate_path("input.csv") == "input.csv"

    def test_traversal_blocked(self):
        with pytest.raises(ValueError, match="Path traversal"):
            _validate_path("../etc/passwd")

    def test_traversal_with_prefix_blocked(self):
        with pytest.raises(ValueError, match="Path traversal"):
            _validate_path("foo/../../../etc/shadow")


class TestProcessBatchErrorHandling:
    """Tests for batch processing error handling."""

    def test_missing_input_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            process_batch(str(tmp_path / "nonexistent.csv"), str(tmp_path / "out.csv"))

    def test_batch_output_content(self, tmp_path):
        csv_in = tmp_path / "in.csv"
        csv_out = tmp_path / "out.csv"
        csv_in.write_text("Patient_ID,v1,v2\nPT-001,10.0,5.0\n", encoding="utf-8")
        process_batch(str(csv_in), str(csv_out))
        content = csv_out.read_text(encoding="utf-8")
        assert "PT-001" in content
        assert "score" in content
        assert "classification" in content


class TestAuditTrailSecurity:
    """Tests for HMAC audit trail security."""

    def test_no_hardcoded_default_key(self):
        """AuditTrail should warn when no key is provided."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            trail = AuditTrail()
            assert len(w) == 1
            assert "AUDIT_SECRET_KEY not set" in str(w[0].message)

    def test_explicit_key_no_warning(self):
        """AuditTrail with explicit key should not warn."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            trail = AuditTrail(secret_key="test-key-12345")
            assert len(w) == 0

    def test_env_var_key_no_warning(self):
        """AuditTrail with env var key should not warn."""
        os.environ["AUDIT_SECRET_KEY"] = "env-test-key-12345"
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                trail = AuditTrail()
                assert len(w) == 0
        finally:
            del os.environ["AUDIT_SECRET_KEY"]

    def test_audit_integrity_with_key(self):
        """Audit trail should maintain integrity with a proper key."""
        trail = AuditTrail(secret_key="test-key-for-integrity")
        trail.log("test", "tier1", "EVENT", {"data": "value1"})
        trail.log("test", "tier1", "EVENT", {"data": "value2"})
        assert trail.verify_integrity() is True
        assert len(trail.get_trail()) == 2


class TestPHIGuard:
    """Tests for PHI outbound guard."""

    def test_clean_text_passes(self):
        PHIGuard.assert_no_phi("Analytical assay specimen KEY-001 optimal")

    def test_mrn_blocked(self):
        with pytest.raises(SecurityException):
            PHIGuard.assert_no_phi("Patient MRN-12345678")

    def test_ssn_blocked(self):
        with pytest.raises(SecurityException):
            PHIGuard.assert_no_phi("SSN: 123-45-6789")

    def test_email_blocked(self):
        with pytest.raises(SecurityException):
            PHIGuard.assert_no_phi("Contact patient@hospital.com")

    def test_redaction(self):
        result = PHIGuard.redact_phi("Patient MRN-12345678 is here")
        assert "REDACTED_IDENTIFIER" in result
        assert "MRN-12345678" not in result


class TestCalculateMetrics:
    """Tests for core calculation function."""

    def test_basic_calculation(self):
        res = calculate_metrics(v1=10.0, v2=5.0)
        assert res["score"] == 12.5  # 10 + 5*(1/2) = 12.5
        assert res["classification"] == "Moderate / Intermediate"

    def test_low_tier(self):
        res = calculate_metrics(v1=5.0)
        assert res["classification"] == "Low / Standard"

    def test_high_tier(self):
        res = calculate_metrics(v1=30.0)
        assert res["classification"] == "High / Severe"

    def test_empty_input(self):
        res = calculate_metrics()
        assert res["score"] == 1.0  # default primary_val
        assert res["inputs_evaluated"] == 0

    def test_string_input_ignored(self):
        res = calculate_metrics(v1=10.0, name="test")
        assert res["score"] == 10.0
        assert res["inputs_evaluated"] == 2
