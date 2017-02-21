import numpy as np
import matplotlib.pyplot as plt

binned_x = np.linspace(-np.pi, np.pi, 101)[:-1]
CA3_x_offset = np.linspace(-np.pi, np.pi, 201)[:-1]
CA3_x_offset_T = CA3_x_offset[np.newaxis].T

field_amplitude = 2.
alpha = 2. # regularization parameter


CA3_rate = (np.exp(field_amplitude * np.cos(binned_x - CA3_x_offset_T)) - np.exp(-field_amplitude)) / \
            (np.exp(field_amplitude) - np.exp(-field_amplitude))
vm = np.exp(field_amplitude * np.cos(binned_x)) + np.exp(field_amplitude * np.cos(binned_x+1.2)) + \
     np.random.uniform(0,1,binned_x.shape)
# vm -= np.min(vm)

# invert CA3_rate matrix
[U,s,Vh] = np.linalg.svd(CA3_rate)
V = Vh.T
D = np.zeros_like(CA3_rate)

D[np.where(np.eye(*D.shape))] = s / (s**2. + alpha**2.)
CA3_inv = V.dot(D.conj().T).dot(U.conj().T)

weights = vm.dot(CA3_inv)
model_vm = weights.dot(CA3_rate)

baseline_indexes = np.where(vm <= np.percentile(vm, 10.))[0]
vm_baseline = np.mean(vm[baseline_indexes])

model_baseline_indexes = np.where(model_vm <= np.percentile(model_vm, 10.))[0]
model_vm_baseline = np.mean(model_vm[model_baseline_indexes])

plt.plot(binned_x, np.subtract(vm, vm_baseline))
plt.plot(binned_x, np.subtract(model_vm, model_vm_baseline))
plt.show()
plt.close()