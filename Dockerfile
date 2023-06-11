# Start with OpenSUSE Tumbleweed
FROM opensuse/tumbleweed

# Install Python 3.10 and other necessary packages
RUN zypper ref && \
    zypper -n in python310 python310-pip zsh uuid-runtime shadow git git-lfs && \
    zypper clean

WORKDIR /app

# Create a non-root user
RUN useradd -ms /bin/bash nyx-user

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

# Install workflow and pipeline features
RUN pip install kedro kedro-viz ruff jupyter

# Make this available
ENV PATH="${PATH}:/home/nyx-user/.local/bin"

# Copy over the wheel files
COPY dist/*.whl /app

# Install the Python wheel
RUN pip install /app/*.whl


CMD ["zsh"]
