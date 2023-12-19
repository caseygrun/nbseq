import pandas as pd
import numpy as np
import sklearn
import skbio.diversity
from skbio import OrdinationResults

from sklearn.decomposition import TruncatedSVD
from sklearn.manifold import MDS
from sklearn.manifold import TSNE


def _ordinate_decomposition_tsne(decomposition_cls = TruncatedSVD, tsne_cls=TSNE):
    def _ordinate(data, samples, observations, n_components=2, decomp_kwargs={}, n_jobs=1, **kwargs):
        method_name = (decomposition_cls.__name__ + '-' + tsne_cls.__name__)

        _ord = decomposition_cls(**decomp_kwargs)
        _ord.fit(data)
        ft_ord = _ord.transform(data)

        _tsne = tsne_cls(n_components=n_components, n_jobs=n_jobs, **kwargs)
        _tsne.fit(ft_ord)
        ft_tsne = _tsne.fit_transform(ft_ord)

        ord_res = OrdinationResults(
            short_method_name=method_name,
            long_method_name=method_name,
            eigvals=pd.Series(np.ones(n_components)/n_components),
            samples=pd.DataFrame(ft_tsne, index=samples),
            proportion_explained=pd.Series(np.ones(n_components)/n_components)
        )
        return _tsne, ord_res
    return _ordinate

def _ordinate_decomposition(method_name, method_cls):
    def _ordinate(data, samples, observations, n_jobs=1, **kwargs):
        _ord = method_cls(**kwargs)
        _ord.fit(data)
        ft_ord = _ord.transform(data)

        ord_res = OrdinationResults(
            short_method_name=method_name,
            long_method_name=_ord.__class__.__name__,
            eigvals=pd.Series(_ord.explained_variance_),
            samples=pd.DataFrame(ft_ord, index=samples),
            features=pd.DataFrame(_ord.components_.T, index=observations),
            proportion_explained=pd.Series(_ord.explained_variance_ratio_)
        )

        return (_ord, ord_res)
    return _ordinate


def mds_eigs(A):
    """https://stackoverflow.com/a/40534620/4091874"""

    # square it
    A = A**2

    # centering matrix
    n = A.shape[0]
    J_c = 1./n*(np.eye(n) - 1 + (n-1)*np.eye(n))

    # perform double centering
    B = -0.5*(J_c.dot(A)).dot(J_c)

    # find eigenvalues and eigenvectors
    eigen_val = la.eig(B)[0]
    #     eigen_vec = la.eig(B)[1].T

    # select top 2 dimensions (for example)
    #     PC1 = np.sqrt(eigen_val[0])*eigen_vec[0]
    #     PC2 = np.sqrt(eigen_val[1])*eigen_vec[1]

    return np.sqrt(eigen_val)

def _ordinate_mds(method_name, method_cls):
    def _ordinate(data, samples, observations, n_jobs=1, **kwargs):
        _ord = method_cls(n_jobs=n_jobs, **kwargs)
        _ord.fit(data)
        ft_ord = _ord.transform(data)

        eigs = mds_eigs(_ord.dissimilarity_matrix_)
        proportion_explained = eigs / eigs.sum()
        ord_res = OrdinationResults(
            short_method_name=method_name,
            long_method_name=_ord.__class__.__name__,
            eigvals=pd.Series(eigs),
            samples=pd.DataFrame(ft_ord, index=samples),
            proportion_explained=pd.Series(proportion_explained)
        )

        return (_ord, ord_res)
    return _ordinate


methods = {
    'tsvd':_ordinate_decomposition('TSVD',TruncatedSVD),
    'mds': _ordinate_mds('MDS',MDS),
    'tsvd-tsne':_ordinate_decomposition_tsne()
}


def ordinate(ft, method, **kwargs):
    import biom
    _ord_f = methods[method.lower()]

    if isinstance(ft, biom.Table):
        return _ord_f(ft.matrix_data.T, samples=ft.ids('sample'), observations=ft.ids('observation'), **kwargs)
    else:
        from anndata import AnnData
        if isinstance(ft, AnnData):
            return _ord_f(ft.X,
                samples=ft.obs_names.values,
                observations=ft.var_names.values, **kwargs)
