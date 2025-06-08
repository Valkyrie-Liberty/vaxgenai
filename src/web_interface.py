"""
Web Interface Module for VaxGenAI

This module provides a web interface for the VaxGenAI platform using FastAPI.
"""

import os
import logging
import tempfile
import shutil
from pathlib import Path
from typing import List, Optional

from fastapi import FastAPI, File, UploadFile, Form, HTTPException
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel

from .utils import get_data_path, get_results_path

logger = logging.getLogger("vaxgenai.web_interface")

class WebInterface:
    """
    Web interface for the VaxGenAI platform
    """
    
    def __init__(self):
        """Initialize the web interface"""
        self.app = FastAPI(title="VaxGenAI", description="AI-Powered Vaccine Design Platform")
        self.data_path = get_data_path()
        self.results_path = get_results_path()
        
        # Create static directory
        self.static_path = Path(__file__).parent.parent / "static"
        self.static_path.mkdir(exist_ok=True)
        
        # Create templates directory
        self.templates_path = Path(__file__).parent.parent / "templates"
        self.templates_path.mkdir(exist_ok=True)
        
        # Create templates
        self.create_templates()
        
        # Mount static files
        self.app.mount("/static", StaticFiles(directory=str(self.static_path)), name="static")
        
        # Set up templates
        self.templates = Jinja2Templates(directory=str(self.templates_path))
        
        # Set up routes
        self.setup_routes()
        
        logger.info("Web interface initialized")
    
    def create_templates(self):
        """Create HTML templates for the web interface"""
        # Create index.html
        index_html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>VaxGenAI - AI-Powered Vaccine Design Platform</title>
            <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css">
            <style>
                body { padding-top: 20px; }
                .container { max-width: 800px; }
                .header { margin-bottom: 30px; text-align: center; }
                .form-group { margin-bottom: 15px; }
                .results { margin-top: 30px; }
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h1>VaxGenAI</h1>
                    <p class="lead">AI-Powered Vaccine Design Platform for Emerging Pathogens</p>
                </div>
                
                <div class="card">
                    <div class="card-header">
                        <h3>Upload Pathogen Sequence</h3>
                    </div>
                    <div class="card-body">
                        <form action="/analyze" method="post" enctype="multipart/form-data">
                            <div class="form-group">
                                <label for="sequence_file">Sequence File (FASTA format):</label>
                                <input type="file" class="form-control" id="sequence_file" name="sequence_file" required>
                            </div>
                            
                            <div class="form-group">
                                <label for="sequence_type">Sequence Type:</label>
                                <select class="form-control" id="sequence_type" name="sequence_type">
                                    <option value="protein">Protein</option>
                                    <option value="dna">DNA</option>
                                    <option value="rna">RNA</option>
                                </select>
                            </div>
                            
                            <div class="form-group">
                                <label for="vaccine_types">Vaccine Types to Design:</label>
                                <div class="form-check">
                                    <input class="form-check-input" type="checkbox" id="subunit" name="vaccine_types" value="subunit" checked>
                                    <label class="form-check-label" for="subunit">Subunit Vaccine</label>
                                </div>
                                <div class="form-check">
                                    <input class="form-check-input" type="checkbox" id="mrna" name="vaccine_types" value="mrna" checked>
                                    <label class="form-check-label" for="mrna">mRNA Vaccine</label>
                                </div>
                                <div class="form-check">
                                    <input class="form-check-input" type="checkbox" id="peptide" name="vaccine_types" value="peptide" checked>
                                    <label class="form-check-label" for="peptide">Peptide Vaccine</label>
                                </div>
                            </div>
                            
                            <button type="submit" class="btn btn-primary">Analyze and Design Vaccines</button>
                        </form>
                    </div>
                </div>
                
                <div class="card mt-4">
                    <div class="card-header">
                        <h3>Example Sequences</h3>
                    </div>
                    <div class="card-body">
                        <p>You can use one of these example sequences:</p>
                        <ul>
                            <li><a href="/static/examples/spike_protein.fasta">SARS-CoV-2 Spike Protein</a></li>
                            <li><a href="/static/examples/influenza_ha.fasta">Influenza Hemagglutinin</a></li>
                        </ul>
                    </div>
                </div>
            </div>
            
            <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
        </body>
        </html>
        """
        
        # Create results.html
        results_html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>VaxGenAI - Results</title>
            <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css">
            <style>
                body { padding-top: 20px; }
                .container { max-width: 800px; }
                .header { margin-bottom: 30px; text-align: center; }
                .results { margin-top: 30px; }
                .nav-tabs { margin-bottom: 15px; }
                .plot-container { text-align: center; margin: 20px 0; }
                .plot-container img { max-width: 100%; }
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h1>VaxGenAI Results</h1>
                    <p class="lead">Analysis and Vaccine Design Results</p>
                </div>
                
                <div class="card">
                    <div class="card-header">
                        <h3>Summary</h3>
                    </div>
                    <div class="card-body">
                        <p><strong>Sequence:</strong> {{ sequence_name }}</p>
                        <p><strong>Sequence Type:</strong> {{ sequence_type }}</p>
                        <p><strong>Sequence Length:</strong> {{ sequence_length }} {{ "amino acids" if sequence_type == "protein" else "nucleotides" }}</p>
                        <p><strong>Epitopes Predicted:</strong> {{ num_epitopes }}</p>
                        <p><strong>Vaccine Candidates Designed:</strong> {{ num_candidates }}</p>
                    </div>
                </div>
                
                <ul class="nav nav-tabs mt-4" id="resultTabs" role="tablist">
                    <li class="nav-item" role="presentation">
                        <button class="nav-link active" id="epitopes-tab" data-bs-toggle="tab" data-bs-target="#epitopes" type="button" role="tab" aria-controls="epitopes" aria-selected="true">Epitopes</button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link" id="vaccines-tab" data-bs-toggle="tab" data-bs-target="#vaccines" type="button" role="tab" aria-controls="vaccines" aria-selected="false">Vaccine Candidates</button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link" id="visualizations-tab" data-bs-toggle="tab" data-bs-target="#visualizations" type="button" role="tab" aria-controls="visualizations" aria-selected="false">Visualizations</button>
                    </li>
                </ul>
                
                <div class="tab-content" id="resultTabsContent">
                    <div class="tab-pane fade show active" id="epitopes" role="tabpanel" aria-labelledby="epitopes-tab">
                        <h4>Predicted Epitopes</h4>
                        <p>Top {{ top_epitopes|length }} epitopes predicted from the sequence:</p>
                        
                        <table class="table table-striped">
                            <thead>
                                <tr>
                                    <th>Position</th>
                                    <th>Sequence</th>
                                    <th>Type</th>
                                    <th>Score</th>
                                </tr>
                            </thead>
                            <tbody>
                                {% for epitope in top_epitopes %}
                                <tr>
                                    <td>{{ epitope.start }}-{{ epitope.end }}</td>
                                    <td>{{ epitope.sequence }}</td>
                                    <td>{{ epitope.type }}</td>
                                    <td>{{ "%.3f"|format(epitope.score) }}</td>
                                </tr>
                                {% endfor %}
                            </tbody>
                        </table>
                        
                        <div class="d-grid gap-2">
                            <a href="/download/epitopes" class="btn btn-primary">Download All Epitopes</a>
                        </div>
                    </div>
                    
                    <div class="tab-pane fade" id="vaccines" role="tabpanel" aria-labelledby="vaccines-tab">
                        <h4>Vaccine Candidates</h4>
                        <p>Top vaccine candidates ranked by overall score:</p>
                        
                        <table class="table table-striped">
                            <thead>
                                <tr>
                                    <th>Rank</th>
                                    <th>Type</th>
                                    <th>Efficacy</th>
                                    <th>Population Coverage</th>
                                    <th>Manufacturability</th>
                                    <th>Overall Score</th>
                                </tr>
                            </thead>
                            <tbody>
                                {% for i, candidate in enumerate(ranked_candidates) %}
                                <tr>
                                    <td>{{ i+1 }}</td>
                                    <td>{{ candidate.type }}</td>
                                    <td>{{ "%.3f"|format(candidate.efficacy_score) }}</td>
                                    <td>{{ "%.3f"|format(candidate.population_coverage) }}</td>
                                    <td>{{ "%.3f"|format(candidate.manufacturability) }}</td>
                                    <td>{{ "%.3f"|format(candidate.overall_score) }}</td>
                                </tr>
                                {% endfor %}
                            </tbody>
                        </table>
                        
                        <div class="d-grid gap-2">
                            <a href="/download/vaccines" class="btn btn-primary">Download Vaccine Candidates</a>
                            <a href="/download/report" class="btn btn-secondary">Download Full Report</a>
                        </div>
                    </div>
                    
                    <div class="tab-pane fade" id="visualizations" role="tabpanel" aria-labelledby="visualizations-tab">
                        <h4>Visualizations</h4>
                        
                        {% for plot in plots %}
                        <div class="plot-container">
                            <h5>{{ plot.title }}</h5>
                            <img src="{{ plot.path }}" alt="{{ plot.title }}">
                        </div>
                        {% endfor %}
                        
                        <div class="d-grid gap-2">
                            <a href="/download/visualizations" class="btn btn-primary">Download All Visualizations</a>
                        </div>
                    </div>
                </div>
                
                <div class="d-grid gap-2 mt-4">
                    <a href="/" class="btn btn-secondary">Start New Analysis</a>
                </div>
            </div>
            
            <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
        </body>
        </html>
        """
        
        # Save templates
        with open(self.templates_path / "index.html", "w") as f:
            f.write(index_html)
        
        with open(self.templates_path / "results.html", "w") as f:
            f.write(results_html)
        
        # Create examples directory
        examples_path = self.static_path / "examples"
        examples_path.mkdir(exist_ok=True)
        
        # Copy spike protein sequence to examples
        shutil.copy(self.data_path / "spike_protein.fasta", examples_path / "spike_protein.fasta")
        
        logger.info("Created web interface templates")
    
    def setup_routes(self):
        """Set up routes for the web interface"""
        @self.app.get("/", response_class=HTMLResponse)
        async def index():
            """Render the index page"""
            return self.templates.TemplateResponse("index.html", {"request": {}})
        
        @self.app.post("/analyze")
        async def analyze(
            sequence_file: UploadFile = File(...),
            sequence_type: str = Form(...),
            vaccine_types: List[str] = Form(...)
        ):
            """
            Analyze a sequence and design vaccines
            
            Args:
                sequence_file: Uploaded sequence file
                sequence_type: Type of sequence (protein, dna, rna)
                vaccine_types: Types of vaccines to design
                
            Returns:
                dict: Analysis results
            """
            try:
                # Save uploaded file
                file_path = self.data_path / sequence_file.filename
                with open(file_path, "wb") as f:
                    f.write(await sequence_file.read())
                
                # TODO: Implement actual analysis and vaccine design
                
                # For now, return a placeholder response
                return {
                    "status": "success",
                    "message": "Analysis completed successfully",
                    "sequence_file": sequence_file.filename,
                    "sequence_type": sequence_type,
                    "vaccine_types": vaccine_types,
                    "results_url": "/results/12345"  # Placeholder
                }
            
            except Exception as e:
                logger.error(f"Error analyzing sequence: {e}")
                raise HTTPException(status_code=500, detail=str(e))
        
        @self.app.get("/results/{job_id}", response_class=HTMLResponse)
        async def results(job_id: str):
            """
            Render the results page
            
            Args:
                job_id: Job ID
                
            Returns:
                HTMLResponse: Results page
            """
            # TODO: Implement actual results retrieval
            
            # For now, return a placeholder response
            return self.templates.TemplateResponse(
                "results.html",
                {
                    "request": {},
                    "job_id": job_id,
                    "sequence_name": "Example Sequence",
                    "sequence_type": "protein",
                    "sequence_length": 1273,
                    "num_epitopes": 42,
                    "num_candidates": 3,
                    "top_epitopes": [
                        {"start": 10, "end": 25, "sequence": "EXAMPLE_EPITOPE_1", "type": "B-cell", "score": 0.95},
                        {"start": 100, "end": 115, "sequence": "EXAMPLE_EPITOPE_2", "type": "T-cell", "score": 0.92},
                        {"start": 200, "end": 215, "sequence": "EXAMPLE_EPITOPE_3", "type": "B-cell", "score": 0.89}
                    ],
                    "ranked_candidates": [
                        {"type": "subunit", "efficacy_score": 0.92, "population_coverage": 0.85, "manufacturability": 0.78, "overall_score": 0.88},
                        {"type": "mRNA", "efficacy_score": 0.88, "population_coverage": 0.82, "manufacturability": 0.75, "overall_score": 0.84},
                        {"type": "peptide", "efficacy_score": 0.85, "population_coverage": 0.80, "manufacturability": 0.82, "overall_score": 0.83}
                    ],
                    "plots": [
                        {"title": "Epitope Scores", "path": "/static/examples/epitope_scores.png"},
                        {"title": "Vaccine Candidate Comparison", "path": "/static/examples/vaccine_comparison.png"}
                    ]
                }
            )
        
        @self.app.get("/download/{file_type}")
        async def download(file_type: str):
            """
            Download files
            
            Args:
                file_type: Type of file to download
                
            Returns:
                FileResponse: File to download
            """
            # TODO: Implement actual file downloads
            
            # For now, return a placeholder response
            return {"status": "success", "message": f"Download {file_type} not yet implemented"}
    
    def run(self, host="0.0.0.0", port=8000):
        """
        Run the web interface
        
        Args:
            host: Host to run on
            port: Port to run on
        """
        import uvicorn
        logger.info(f"Starting web interface on {host}:{port}")
        uvicorn.run(self.app, host=host, port=port)

