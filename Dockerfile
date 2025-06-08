# Use Python 3.9 slim image as base
FROM python:3.9-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt requirements-dev.txt ./

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY src/ ./src/
COPY data/ ./data/
COPY models/ ./models/
COPY static/ ./static/
COPY templates/ ./templates/
COPY *.py ./

# Create necessary directories
RUN mkdir -p results/visualizations results/reports

# Set Python path
ENV PYTHONPATH=/app

# Expose port for web interface
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Default command
CMD ["python", "run_web.py", "--host", "0.0.0.0", "--port", "8000"]

# Alternative commands for different use cases:
# For demo: CMD ["python", "demo.py"]
# For testing: CMD ["python", "comprehensive_testing.py"]
# For batch processing: CMD ["python", "-m", "src.batch_processor"]

