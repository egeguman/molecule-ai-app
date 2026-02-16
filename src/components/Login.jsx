import React, { useState } from 'react';

const Login = ({ onLogin, onSwitchToRegister }) => {
    const [formData, setFormData] = useState({
        email: '',
        password: ''
    });
    const [error, setError] = useState('');
    const [loading, setLoading] = useState(false);

    const handleChange = (e) => {
        setFormData({
            ...formData,
            [e.target.name]: e.target.value
        });
        if (error) setError('');
    };

    const handleSubmit = async (e) => {
        e.preventDefault();
        setError('');
        setLoading(true);

        try {
            const response = await fetch('/api/auth/login', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    email: formData.email,
                    password: formData.password,
                }),
            });
            const data = await response.json();

            if (!response.ok) {
                setError(data.error || 'Login failed. Please try again.');
                return;
            }

            localStorage.setItem('token', data.token);
            if (onLogin) onLogin(data.user);
        } catch (err) {
            setError('Unable to connect to server. Please try again.');
        } finally {
            setLoading(false);
        }
    };

    return (
        <div style={styles.container}>
            <div style={styles.card}>
                <div style={styles.header}>
                    <h2 style={styles.title}>Welcome Back</h2>
                    <p style={styles.subtitle}>Sign in to continue your research</p>
                </div>

                <form onSubmit={handleSubmit} style={styles.form}>
                    {error && (
                        <div style={styles.error}>{error}</div>
                    )}

                    <div style={styles.inputGroup}>
                        <label style={styles.label}>Email Address</label>
                        <input
                            type="email"
                            name="email"
                            value={formData.email}
                            onChange={handleChange}
                            placeholder="Enter your email"
                            style={styles.input}
                            required
                            disabled={loading}
                        />
                    </div>

                    <div style={styles.inputGroup}>
                        <label style={styles.label}>Password</label>
                        <input
                            type="password"
                            name="password"
                            value={formData.password}
                            onChange={handleChange}
                            placeholder="Enter your password"
                            style={styles.input}
                            required
                            disabled={loading}
                        />
                    </div>

                    <button type="submit" style={{
                        ...styles.button,
                        opacity: loading ? 0.7 : 1,
                        cursor: loading ? 'not-allowed' : 'pointer',
                    }} disabled={loading}>
                        {loading ? 'Signing in...' : 'Sign In'}
                    </button>
                </form>

                <div style={styles.switchSection}>
                    <span style={styles.switchText}>Don't have an account?</span>
                    <button
                        type="button"
                        style={styles.switchButton}
                        onClick={onSwitchToRegister}
                    >
                        Create Account
                    </button>
                </div>
            </div>
        </div>
    );
};

const styles = {
    container: {
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        minHeight: '600px',
        padding: '20px',
        animation: 'fadeIn 0.5s ease-out',
    },
    card: {
        width: '100%',
        maxWidth: '400px',
        backgroundColor: '#ffffff',
        borderRadius: '16px',
        boxShadow: '0 4px 20px rgba(0, 0, 0, 0.05)',
        padding: '40px',
        border: '1px solid #f3f4f6',
    },
    header: {
        textAlign: 'center',
        marginBottom: '32px',
    },
    title: {
        fontSize: '28px',
        fontWeight: '700',
        color: '#111827',
        margin: '0 0 8px 0',
        fontFamily: "'DM Sans', sans-serif",
    },
    subtitle: {
        fontSize: '15px',
        color: '#6b7280',
        margin: 0,
        lineHeight: '1.5',
        fontFamily: "'DM Sans', sans-serif",
    },
    error: {
        padding: '10px 14px',
        backgroundColor: '#fef2f2',
        border: '1px solid #fecaca',
        borderRadius: '8px',
        color: '#dc2626',
        fontSize: '13px',
        fontFamily: "'DM Sans', sans-serif",
    },
    form: {
        display: 'flex',
        flexDirection: 'column',
        gap: '20px',
    },
    inputGroup: {
        display: 'flex',
        flexDirection: 'column',
        gap: '8px',
    },
    label: {
        fontSize: '14px',
        fontWeight: '600',
        color: '#374151',
        fontFamily: "'DM Sans', sans-serif",
    },
    input: {
        padding: '12px 16px',
        borderRadius: '8px',
        border: '1px solid #e5e7eb',
        fontSize: '15px',
        backgroundColor: '#f9fafb',
        transition: 'all 0.2s ease',
        outline: 'none',
        fontFamily: "'DM Sans', sans-serif",
    },
    button: {
        marginTop: '12px',
        padding: '14px',
        backgroundColor: '#16a34a',
        color: '#ffffff',
        border: 'none',
        borderRadius: '8px',
        fontSize: '16px',
        fontWeight: '600',
        cursor: 'pointer',
        transition: 'background-color 0.2s ease',
        fontFamily: "'DM Sans', sans-serif",
    },
    switchSection: {
        display: 'flex',
        justifyContent: 'center',
        alignItems: 'center',
        gap: '8px',
        marginTop: '24px',
        paddingTop: '20px',
        borderTop: '1px solid #f3f4f6',
    },
    switchText: {
        fontSize: '14px',
        color: '#6b7280',
        fontFamily: "'DM Sans', sans-serif",
    },
    switchButton: {
        fontSize: '14px',
        fontWeight: '600',
        color: '#16a34a',
        background: 'none',
        border: 'none',
        cursor: 'pointer',
        fontFamily: "'DM Sans', sans-serif",
        padding: 0,
        transition: 'color 0.2s ease',
    }
};

export default Login;
