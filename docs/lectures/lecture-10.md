<script type="text/x-thebe-config">
  {
      requestKernel: true,
      mountActivateWidget: true,
      mountStatusWidget: true,
      binderOptions: {
      repo: "mmosmond/executable-cells",
      ref: "main",
      },
  }
</script>
<script src="https://unpkg.com/thebe@latest/lib/index.js"></script>
<link rel="stylesheet" href="https://unpkg.com/thebe@latest/lib/thebe.css">

# Lecture 10: Linear algebra I

<hr style="margin-bottom: 0em;">
<center>
<div class="inrow">
	Run notes interactively?
	<div style="float: left;" class="thebe-activate"></div>
	<div class="thebe-status"></div>
</div>
</center>
<hr style="margin-top: 0em;">

## Lecture overview

1. [Motivation](#section1)
2. [What are vectors?](#section2)
3. [What is a matrix?](#section3)
4. [Vector and matrix operations](#section4)

<span id='section1'></span>
## 1. Motivation
<hr>

Until now, we have been dealing with problems in a single variable changing over time.
  
Often, dynamical systems involve more than one variable (ie, they are **multivariate**). For instance, we may be interested in how the numbers of two species change as they interact (e.g., compete) with one another.

As a simple example with more than one variable, consider a model tracking the number of birds on two islands. Let the number of birds on island 1 be $n_1$ and let the number of birds on island 2 be $n_2$. We assume the birds migrate between the islands at per capita rates $m_{12}$ and $m_{21}$, the birds on each island give birth at per capita rates $b_1$ and $b_2$, the birds on each island die at per capita rates $d_1$ and $d_2$, and new birds arrive on each island at rates $m_1$ and $m_2$. This is captured in the following flow diagram

<center>
```mermaid
graph LR;
    A1((n1)) --b1 n1--> A1;
    B1[ ] --m1--> A1;
    A1 --d1 n1--> C1[ ];

    A2((n2)) --b2 n2--> A2;
    B2[ ] --m2--> A2;
    A2 --d2 n2--> C2[ ];

    A1 --m12 n1--> A2;
    A2 --m21 n2--> A1;

    style B1 height:0px;
    style C1 height:0px;
    style B2 height:0px;
    style C2 height:0px;
```   
</center>
    
The rate of change in $n_1$ and $n_2$ are then described by the following system of differential equations

$$
\begin{aligned}
\frac{\mathrm{d}n_1}{\mathrm{d}t} &= (b_1 - d_1 - m_{12})n_1 + m_{21} n_2 + m_1 \\
\frac{\mathrm{d}n_1}{\mathrm{d}t} &= m_{12} n_1 + (b_2 - d_2 - m_{21})n_2 + m_1
\end{aligned}
$$

These equations are linear functions of the variables (i.e., they contain only constant multiples of $n_1$ and $n_2$ and nothing more complicated such as $x^2$ or $e^x$).

Linear systems of equations like these can also be written in **matrix form**

$$
\begin{aligned}
\begin{pmatrix} \frac{\mathrm{d}n_1}{\mathrm{d}t} \\ \frac{\mathrm{d}n_2}{\mathrm{d}t} \end{pmatrix} 
&= \begin{pmatrix} b_1 - d_1 - m_{12} & m_{21} \\ m_{12} & b_2 - d_2 - m_{21} \end{pmatrix}
\begin{pmatrix} n_1 \\ n_2 \end{pmatrix} 
+ \begin{pmatrix} m_1 \\ m_2 \end{pmatrix}\\
\frac{\mathrm{d}\vec{n}}{\mathrm{d}t} &= \mathbf{M}\vec{n} + \vec{m}
\end{aligned}
$$

Not only is this a nice compact expression, there are rules of linear algebra that can help us conveniently solve this (and any other) set of linear equations.

So let's get to know these rules.

<span id='section2'></span>
## 2. What are vectors?
<hr>

**Vectors** are lists of elements (elements being numbers, parameters, functions, or variables).

A **column vector** has elements arranged from top to bottom

$$
\begin{equation*}
\begin{pmatrix}
  5 \\
  2
\end{pmatrix},
\begin{pmatrix}
  1 \\
  5 \\
  9 \\
  7
\end{pmatrix},
\begin{pmatrix}
  x \\
  y
\end{pmatrix},
\begin{pmatrix}
  x \\
  y \\
  z
\end{pmatrix},
\begin{pmatrix}
  x_1 \\
  x_2 \\
  \vdots \\
  x_n
\end{pmatrix}
\end{equation*}
$$

A **row vector** has elements arranged from left to right

$$
\begin{pmatrix}5 & 2\end{pmatrix}, \begin{pmatrix} 1 & 5 & 9 & 7\end{pmatrix}, \begin{pmatrix} x & y \end{pmatrix}, \begin{pmatrix} x & y & z\end{pmatrix}, \begin{pmatrix} x_1 & x_2 & \cdots & x_n \end{pmatrix}
$$

We will indicate vectors by placing an arrow on top of the symbol

$$
\vec{x} = \begin{pmatrix} x_1 & x_2 & \cdots & x_n \end{pmatrix}
$$

The number of elements in the vector indicates its **dimension**, $n$.

For example, the row vector $\begin{pmatrix}x & y\end{pmatrix}$ has dimension $n=2$.

You can represent a vector as an arrow in $n$ dimensions, connecting the origin
with a point whose coordinates are given by elements in the vector. For example, the vector $\vec{v} = \begin{pmatrix} 1\\2\end{pmatrix}$ can be depicted as below


<pre data-executable="true" data-language="python">
import matplotlib.pyplot as plt #import plotting library

plt.arrow(0, 0, #starting x and y values of arrow
          1, 2, #change in x and y 
          head_width=0.1, color='black') #aesthetics
plt.xlim(0,2.5) #set bounds on x axis
plt.ylim(0,2.5) #set bounds on y axis
plt.show()
</pre>


    
![png](lecture-10_files/lecture-10_3_0.png)
    


<span id='section3'></span>
## 3. What is a matrix?
<hr>

An $m \times n$ **matrix** has $m$ rows and $n$ columns

$$
\begin{equation*}
\begin{pmatrix}
  x_{11}  & x_{12} & \cdots & x_{1n}\\
  x_{21}  & x_{22} & \cdots & x_{2n}\\
  \vdots & \vdots &        & \vdots\\
  x_{m1}  & x_{m2} & \cdots & x_{mn}\\
\end{pmatrix},
\begin{pmatrix}
  a & b \\
  c & d
\end{pmatrix},
\begin{pmatrix}
  75 & 67 \\
  66 & 34 \\
  12 & 14
\end{pmatrix},
\begin{pmatrix}
  1 & 0 & 0 \\
  0 & 1 & 0 \\
  0 & 0 & 1
\end{pmatrix}
\end{equation*}
$$

We will indicate matrices by bolding the symbol (and using capital letters)

$$
\mathbf{X} = \begin{pmatrix}
  x_{11}  & x_{12} & \cdots & x_{1n}\\
  x_{21}  & x_{22} & \cdots & x_{2n}\\
  \vdots & \vdots &        & \vdots\\
  x_{m1}  & x_{m2} & \cdots & x_{mn}\\
\end{pmatrix}
$$

A matrix with an equal number of rows and columns, $m=n$, is a **square matrix**.
  
A matrix with zeros everywhere except along the diagonal is called a **diagonal matrix**
  
$$
\begin{pmatrix}
  a & 0 & 0 \\
  0 & b & 0 \\
  0 & 0 & c
\end{pmatrix}
$$
   
And a special case of this with 1s along the diagonal is called the **identity matrix**

$$
\begin{pmatrix}
  1 & 0 & 0 \\
  0 & 1 & 0 \\
  0 & 0 & 1
\end{pmatrix}
$$
  
A matrix with all zeros below the diagonal is called an **upper trianglular matrix**

$$
\begin{pmatrix}
  a & b & c \\
  0 & d & e \\
  0 & 0 & f
\end{pmatrix}
$$
  
A matrix with all zeros above the diagonal is called an **lower trianglular matrix**

$$
\begin{pmatrix}
  a & 0 & 0 \\
  b & d & 0 \\
  c & e & f
\end{pmatrix}
$$

It is sometimes useful to chop a matrix up into multiple blocks, creating a **block matrix**

$$
\begin{pmatrix}
  a & b & c \\
  d & e & f \\
  g & h & i
\end{pmatrix}
= 
\begin{pmatrix}
  \mathbf{A} & \mathbf{B} \\
  \mathbf{C} & \mathbf{D}
\end{pmatrix}
$$

where $\mathbf{A}=\begin{pmatrix} a & b \\ d & e\end{pmatrix}$, $\mathbf{B}=\begin{pmatrix} c\\ f\end{pmatrix}$, $\mathbf{C}=\begin{pmatrix} g & h \end{pmatrix}$, and $\mathbf{D}=\begin{pmatrix} i\end{pmatrix}$.

This is especially helpful when the block form has off-diagonal submatrices consisting of all zeros.

For instance, when $\mathbf{B}=\begin{pmatrix} 0\\0\end{pmatrix}$ or $\mathbf{C}=\begin{pmatrix} 0 & 0 \end{pmatrix}$, we have a **block triangular matrix**.

And when $\mathbf{B}=\begin{pmatrix} 0\\0\end{pmatrix}$ and $\mathbf{C}=\begin{pmatrix} 0 & 0 \end{pmatrix}$, 
we have a **block diagonal matrix**.

Finally, it is sometimes useful to **transpose** a matrix, which exchanges the rows and columns (an element in row $i$ column $j$ moves to row $j$ column $i$)

$$
\begin{pmatrix}
  a_1 & a_2 & a_3 \\
  b_1 & b_2 & b_3
\end{pmatrix}^\intercal
= 
\begin{pmatrix}
  a_1 & b_1  \\
  a_2 & b_2  \\
  a_3 & b_3
\end{pmatrix}
$$

Like vectors, matrices have a graphical/geometrical interpretation: they stretch and rotate vectors (as we will see shortly).

<span id='section4'></span>
## 4. Vector and matrix operations
<hr>

### Addition

Vector and matrix **addition** (and subtraction) is straightforward, entry-by-entry:

$$
\begin{equation*}
\begin{pmatrix}
  a \\
  b
\end{pmatrix}
+
\begin{pmatrix}
  c \\
  d
\end{pmatrix}
=
\begin{pmatrix}
  a+c \\
  b+d
\end{pmatrix}
\end{equation*}
$$

$$
\begin{equation*}
\begin{pmatrix}
  a & b \\
  c & d
\end{pmatrix}
+
\begin{pmatrix}
  e & f \\
  g & h
\end{pmatrix}
=
\begin{pmatrix}
  a+e & b+f \\
  c+g & d+h
\end{pmatrix}
\end{equation*}
$$

!!! warning 

    The vectors or matrices added together must have the same dimension!
    
Geometrically, adding vectors is like placing the second vector at the end of the first. Below we add the black and red vectors together to get the blue vector.


<pre data-executable="true" data-language="python">
import matplotlib.pyplot as plt #import plotting library

v1 = [1,2] #vector 1
v2 = [1,0] #vector 2
v12 = [i+j for i,j in zip(v1,v2)] #sum of the two vectors

#first vector
plt.arrow(0, 0, #starting x and y values of arrow
          v1[0], v1[1], #change in x and y 
          head_width=0.1, color='black') #aesthetics

#second vector placed at the end of first vector
plt.arrow(v1[0], v1[1], #starting x and y values of arrow
          v2[0], v2[1], #change in x and y 
          head_width=0.1, color='red') #aesthetics

#sum of the vectors
plt.arrow(0, 0, #starting x and y values of arrow
          v12[0], v12[1], #change in x and y 
          head_width=0.1, color='blue') #aesthetics

plt.xlim(0,2.5) #set bounds on x axis
plt.ylim(0,2.5) #set bounds on y axis
plt.show()
</pre>


    
![png](lecture-10_files/lecture-10_6_0.png)
    


### Multiplication

Vector and matrix **multiplication** by a scalar (which may be a constant, a variable, or a function, but not a matrix or a vector) is also straightforward, we just multiply every element by the scalar:

$$
\begin{equation*}
\alpha *
\begin{pmatrix}
  a \\
  b
\end{pmatrix}
=
\begin{pmatrix}
  \alpha a \\
  \alpha b
\end{pmatrix}
\end{equation*}
$$

$$
\begin{equation*}
\alpha *
\begin{pmatrix}
  a & b\\
  c & d
\end{pmatrix}
=
\begin{pmatrix}
  \alpha a & \alpha b\\
  \alpha c & \alpha d
\end{pmatrix}
\end{equation*}
$$

Geometrically, multiplying by a scalar stretches (if $\alpha>1$) or compresses (if $\alpha<1$) a vector. Below we multiply the black vector by $1/2$ to get the red vector.


<pre data-executable="true" data-language="python">
import matplotlib.pyplot as plt #import plotting library

v1 = [1,2] #vector 1
alpha = 1/2 #scalar
v2 = [i*alpha for i in v1] #multiplication by a scalar

#original vector
plt.arrow(0, 0, #starting x and y values of arrow
          v1[0], v1[1], #change in x and y 
          head_width=0.1, color='black') #aesthetics

#stretched vector
plt.arrow(0, 0, #starting x and y values of arrow
          v2[0], v2[1], #change in x and y 
          head_width=0.1, color='red') #aesthetics

plt.xlim(0,2.5) #set bounds on x axis
plt.ylim(0,2.5) #set bounds on y axis
plt.show()
</pre>


    
![png](lecture-10_files/lecture-10_8_0.png)
    


Multiplying vectors and matrices together is a bit trickier, but is based on the fact that a row vector times a column vector is equal to the sum of the products of their respective entries

$$
\begin{equation*}
\begin{pmatrix} a & b & c \end{pmatrix}
\begin{pmatrix}
  x \\
  y \\
  z
\end{pmatrix}
= ax + by + cz
\end{equation*}
$$

This is referred to as the **dot product**. (There are other types of products for vectors and matrices, which we won't cover in this class.)

To multiply a matrix by a vector, this procedure is repeated first for the first row of the matrix, then for the second row of the matrix, etc:

$$
\begin{equation*}
\begin{pmatrix}
  a & b & c \\
  d & e & f \\
  g & h & i
\end{pmatrix}
\begin{pmatrix}
  x \\
  y \\
  z
\end{pmatrix}
=
\begin{pmatrix}
  ax + by + cz \\
  dx + ey + fz \\
  gx + hy + iz \\
\end{pmatrix}
\end{equation*}
$$

Geometrically, multiplying a vector by a matrix stretches and rotates a vector. Below we multiply the black vector my a matrix to get the red vector.


<pre data-executable="true" data-language="python">
import matplotlib.pyplot as plt #import plotting library
from sympy import *

v = Matrix([[2],[1]]) #column vector
M = Matrix([[1,-1],[1,1/4]]) #matrix
u = M*v

#original vector
plt.arrow(0, 0, #starting x and y values of arrow
          float(v[0]), float(v[1]), #change in x and y 
          head_width=0.1, color='black') #aesthetics

#stretched and rotated vector
plt.arrow(0, 0, #starting x and y values of arrow
          float(u[0]), float(u[1]), #change in x and y 
          head_width=0.1, color='red') #aesthetics

plt.xlim(0,2.5) #set bounds on x axis
plt.ylim(0,2.5) #set bounds on y axis
plt.show()
</pre>


    
![png](lecture-10_files/lecture-10_10_0.png)
    


To multiply a matrix by a matrix, this procedure is then repeated first for the first column of the second matrix and then for the second column of the second matrix, etc:

$$
\begin{equation*}
\begin{pmatrix}
  a & b \\
  c & d
\end{pmatrix}
\begin{pmatrix}
  e & f \\
  g & h
\end{pmatrix}
=
\begin{pmatrix}
  ae + bg & af + bh \\
  ce + dg & cf + dh
\end{pmatrix}
\end{equation*}
$$

!!! warning 

    An $m \times n$ matrix (or vector) $\mathbf{A}$ can be multiplied on the right by $\mathbf{B}$ *only* if $\mathbf{B}$ is an $n \times p$ matrix (or vector). The resulting matrix (or vector) will then be $m \times p$.

As opposed to basic algebra, matrix multiplication is *not* commutative. That is, $\mathbf{AB}$ does not generally equal $\mathbf{BA}$.

This means that if we want to multiply both sides of an equation, e.g., $\mathbf{AB} = \mathbf{C}$, by $\mathbf{D}$, we need to do so on the same side, $\mathbf{ABD} = \mathbf{CD}$ or $\mathbf{DAB} = \mathbf{DC}$. We therefore often need to specify that we are multiplying by a matrix "on the left" or "on the right".

On the other hand, matrix multiplication does satisfy the following laws:

- $(\mathbf{AB})\mathbf{C} = \mathbf{A}(\mathbf{BC})$ (associative law)
- $\mathbf{A}(\mathbf{B+C}) = \mathbf{AB}+\mathbf{AC}$ (distributive law)
- $(\mathbf{A}+\mathbf{B})\mathbf{C} = \mathbf{AC}+\mathbf{BC}$ (distributive law)
- $\alpha(\mathbf{AB}) = (\alpha\mathbf{A})\mathbf{B} = \mathbf{A}(\alpha\mathbf{B}) = (\mathbf{A}\mathbf{B})\alpha$ (commutative law for scalars)

Multiplication between the identity matrix and any vector, $\vec{v}$, or square matrix, $\mathbf{M}$, has no effect (it is like a "1" in normal algebra)

$$
\mathbf{I}\vec{v}=\vec{v}
$$

$$
\mathbf{I}\mathbf{M}=\mathbf{M}\mathbf{I}=\mathbf{M}
$$
