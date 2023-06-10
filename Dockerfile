# Dockerfile
FROM python:3.10

WORKDIR /app

# Create a non-root user
RUN useradd -ms /bin/bash nyx-user

# Copy over the wheel files
# COPY dist/*.whl /app
COPY target/wheels/*.whl /app

# Install the Python wheel
RUN pip install /app/*.whl

# Install workflow and pipeline features
RUN pip install kedro kedro-viz ruff jupyter

# Install zsh
RUN apt-get update && apt-get install -y zsh

# Switch to the new user
USER nyx-user

# Install oh-my-zsh
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"

# Set the zsh theme
RUN sed -i 's/ZSH_THEME="robbyrussell"/ZSH_THEME="fino-time"/' ~/.zshrc

# Change default shell to zsh
SHELL ["/bin/zsh", "-c"]

RUN mkdir -p $HOME/workspace

WORKDIR /home/nyx-user/workspace

CMD ["zsh"]
