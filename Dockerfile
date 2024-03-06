# Use an official Python runtime as a parent image
FROM python:3.8

# Set the working directory inside the container
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Clone the OpenPNM repository and install it
RUN git clone https://github.com/PMEAL/OpenPNM \
    && cd OpenPNM \
    && git checkout v2.6.0 \
    && pip install --no-cache-dir -r requirements/pip_requirements.txt \
    && pip install -e .
    
# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

RUN cp patch/* OpenPNM/ -r

# Run the main.py file with parameters
CMD ["python", "GA_main_Windows.py"]
