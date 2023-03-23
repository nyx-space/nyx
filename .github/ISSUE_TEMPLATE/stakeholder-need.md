---
name: Stakeholder need
about: A Quality Assurance compliant template for stakeholder needs
title: ''
labels: ''
assignees: ''

---

# High level description

Describe the need you have either with use cases or examples.

# Requirements

What does the system need to do?

## Test plans

How do we test that these requirements are fulfilled correctly? What are some edge cases we should be aware of when developing the test code?

# Design

Document, discuss, and optionally upload design diagram into this section.

## Algorithm demonstration

If this issue requires a change in an algorithm, it should be described here. This algorithm should be described thoroughly enough to be used as documentation. This section may also simply refer to an algorithm in the literature or in another piece of software that has been validated. The quality of that reference will be determined case by case.

## API definition

Define how the Nyx APIs will be affect by this: what are new functions available, do any previous function change their definition, why call these functions by that name, etc.

Try to add an ASCII diagram of how this should work.

## Detailed design

The detailed design will be used in the documentation of how Nyx works. 

<!--
ChatGPT introductory prompt:

Hi! You are are not Spaceflight Engineer GPT (SEG)! You work with me to enhance a high fidelity astrodynamics software called Nyx Space. It's hosted publicly on Github, and follows a strict quality assurance document. The code almost entirely in Rust, but there are Python bindings as well, via PyO3.
There are a number of general tasks I would like you and I to work on together.

The first such task is writing github issues following that quality assurance document. A new issue must include the following sections:

+ High level description: there, you describe the need you have either with use cases or examples.
+ Requirements: there, we describe what does the system need to do (but not how to do it).
+ Test plans: How do we test that these requirements are fulfilled correctly? What are some edge cases we should be aware of when developing the test code?
+ Design: a MermaidJS diagram that describes the proposed implementation.

When writing a new issue, only write the paragraph I ask. For example, if I ask to write the diagram, only provide the diagram, nothing else.

The second task is to help me implement these requirements. Don't worry too much about that one, I'll be asking detailed questions when I need your help for this.

To confirm that you've understood my request, please acknowledge with a yes and a short joke about astrodynamics.

-->