# VCF Consensus Heterozygosity Rate

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

VCF Transition/Transversion (Ti/Tv) & Heterozygosity Analyzer
Computes Ti/Tv ratio, heterozygous/homozygous variant ratios, and missingness metrics across samples.

Zero-dependency Python implementation with single and batch evaluation.
Author: Dr. Abu Suraih Sakhri
License: MIT

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`calculate_metrics()`**: Core domain algorithm for vcf-consensus-heterozygosity-rate.
- **`process_single()`** — calculates and validates process_single parameters.
- **`process_batch()`** — calculates and validates process_batch parameters.
- **`main()`** — calculates and validates main parameters.

---

## 💻 Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/vcf-consensus-heterozygosity-rate.git
cd vcf-consensus-heterozygosity-rate

# Install dependencies
pip install -e ".[dev]"
```

---

## 💻 CLI Quickstart & Usage

### 1. Single Evaluation
```bash
python vcf_titv.py single --v1 12.0 --v2 4.0 --v3 2.0
```

### 2. Batch Processing
```bash
python vcf_titv.py batch -i sample.csv -o results.csv
```

### 3. Enterprise CLI (Full Feature Set)
```bash
# Audit a single task
python cli.py audit --task-id TASK-001 --primary 28.5 --secondary 14.2

# Batch process CSV records
python cli.py batch -i sample.csv -o results.csv

# Verify HMAC audit trail integrity
python cli.py verify-audit

# Launch FastAPI REST server
python cli.py serve --host 127.0.0.1 --port 8000
```

### Parameter Reference
- `--v1`, `--v2`, `--v3`: Numeric parameters for single evaluation
- `-i`, `--input`: Input CSV file path for batch processing
- `-o`, `--output`: Output CSV file path (default: results.csv)

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Parameter / observation metric | Required |
| `v1` | Parameter / observation metric | Required |
| `v2` | Parameter / observation metric | Required |
| `v3` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

### Security Configuration

Set the `AUDIT_SECRET_KEY` environment variable for persistent audit integrity:

```bash
export AUDIT_SECRET_KEY="your-secure-random-key"
```

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py 1000
```

---

## 🐳 Container Deployment

```bash
docker build -t vcf-consensus-heterozygosity-rate .
docker run -p 8000:8000 vcf-consensus-heterozygosity-rate
```

---

## 📁 Project Structure

```
vcf-consensus-heterozygosity-rate/
├── agents/                 # Enterprise agent framework
│   ├── __init__.py
│   ├── api.py             # FastAPI REST server
│   ├── base.py            # Security, PHI guard, audit trail
│   ├── learning.py        # Bayesian calibration engine
│   ├── llm_factory.py     # LLM provider factory
│   ├── metrics.py         # Prometheus metrics collector
│   ├── models.py          # Pydantic data models
│   ├── streamer.py        # WebSocket telemetry streamer
│   ├── supervisor.py      # Multi-agent orchestrator
│   └── workers.py         # Specialized domain workers
├── tests/                 # Test suite
│   ├── test_enrichment.py
│   └── test_vcf_consensus_heterozygosity_rate.py
├── cli.py                 # Enterprise CLI entry point
├── vcf_titv.py            # Core Ti/Tv analysis module
├── enrichment.py          # Enrichment feature modules
├── simulator.py           # High-throughput simulation
├── sample.csv             # Sample input data
├── pyproject.toml         # Python project configuration
├── Dockerfile             # Container build config
└── docker-compose.yml     # Multi-container orchestration
```
